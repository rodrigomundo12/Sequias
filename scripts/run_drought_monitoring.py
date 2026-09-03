# -*- coding: utf-8 -*-

"""
EL SALVADOR DROUGHT MONITORING
Sentinel-2 L2A

Workflow:

    Sentinel-2 monthly data
            ↓
    monthly spectral indices
            ↓
    quarterly mean
            ↓
    leave-one-year-out anomaly
            ↓
    2-98 percentile filtering
            ↓
    5 anomaly classes
            ↓
    direct polygonization
            ↓
    web-optimized GeoJSON

Quarter definitions:

    Q1 = May, June, July
    Q2 = August, September, October
    Q3 = November, December, January
    Q4 = February, March, April

Indices:

    NDDI
    MNDWI
    NMDI
    NDII
    NDVI

The workflow automatically:

    - determines the latest complete quarter
    - uses the latest 3 hydrological years
    - skips Sentinel Hub downloads when composites already exist
    - masks products to the actual El Salvador polygon
    - calculates leave-one-year-out anomalies
    - removes values outside P2-P98
    - polygonizes directly from anomaly values
    - creates 5 classes
    - creates web GeoJSON products
    - creates metadata.json

Designed to run using the ArcGIS Pro Python environment.
"""


# ================================================================
# 1. IMPORTS
# ================================================================

import os
import math
import calendar
import gc
import json
import re
import gzip
import shutil

from datetime import datetime

import requests
import numpy as np
import geopandas as gpd
import rasterio
import fiona

from tqdm import tqdm

from rasterio.merge import merge
from rasterio.transform import from_bounds, Affine
from rasterio.features import shapes
from rasterio.mask import mask
from affine import Affine
from shapely.geometry import shape, mapping

from sentinelhub import (
    SHConfig,
    SentinelHubRequest,
    MimeType,
    CRS,
    BBox,
    DataCollection,
    BBoxSplitter
)


# ================================================================
# 2. MAIN CONFIGURATION
# ================================================================

BASE_DIR = (
    r"C:\Users\LLRL-PAHs2\Downloads\sequias-poligonos"
)

DATA_DIR = os.path.join(
    BASE_DIR,
    "data"
)

DISTRIBUTION_DIR = os.path.join(
    BASE_DIR,
    "distribution"
)

METADATA_DIR = os.path.join(
    BASE_DIR,
    "metadata"
)

MONTHLY_DIR = os.path.join(
    DATA_DIR,
    "monthly_tiles"
)

QUARTERLY_TILE_DIR = os.path.join(
    DATA_DIR,
    "quarterly_tiles"
)

QUARTERLY_DIR = os.path.join(
    DATA_DIR,
    "quarterly_composites"
)

ANOMALY_DIR = os.path.join(
    DATA_DIR,
    "anomalies"
)

POLYGON_DIR = os.path.join(
    DATA_DIR,
    "polygonized"
)

LATEST_RASTER_DIR = os.path.join(
    DISTRIBUTION_DIR,
    "latest",
    "rasters"
)

LATEST_VECTOR_DIR = os.path.join(
    DISTRIBUTION_DIR,
    "latest",
    "vectors"
)

ARCHIVE_DIR = os.path.join(
    DISTRIBUTION_DIR,
    "archive"
)

CROPPED_QUARTERLY_DIR = os.path.join(
    QUARTERLY_DIR,
    "aoi"
)


for directory in [
    DATA_DIR,
    MONTHLY_DIR,
    QUARTERLY_TILE_DIR,
    QUARTERLY_DIR,
    ANOMALY_DIR,
    POLYGON_DIR,
    LATEST_RASTER_DIR,
    LATEST_VECTOR_DIR,
    ARCHIVE_DIR,
    METADATA_DIR,
    CROPPED_QUARTERLY_DIR,
]:

    os.makedirs(
        directory,
        exist_ok=True
    )


# ================================================================
# 3. SENTINEL HUB CREDENTIALS
# ================================================================
#
# PUT YOUR EXISTING CREDENTIAL TUPLES HERE.
#
# Example:
#
CREDENTIALS = [
    ("6fad6783-033a-4f48-ad65-25b7c2d530f6", "GYI9ceecoGvy5YDVEHOb1m3cVLzT1Ieg"),
    ("4bf7191e-bd75-47be-8c33-6068cd706be0", "hgx4XtR66imiKwWIA493MglvsYi83b5i"),
    ("f618b4e5-7fa2-4b61-86a1-5ecc4b469ffd", "ekjzotsnAFZesT7WtN06tp1aBt9pK3fo"),
    ("4b6356f4-cddc-4c3b-a92a-06689245c899", "YUK7xNbPyAq9g4tJ5GYUtsf7gQp7qgez"),
    ("0f1b50bb-bd7a-4273-bc98-6e6ed8773938", "CR6PjgWH7WVxuI4KKahxgWMAPEg3d2qi"),
    ("b37250a1-9665-45c6-8c4a-15c2d3b6e244", "vDbTctISCwmfk4RjIZV4fLbeKDyB6G"),
    ("2da61b8f-5b6c-40f3-b9da-830069339d5f", "E6Runw6PjHPNOBnoTLRvdcblWI6U5zop"),
    ("4da6eb44-10cf-42f4-9716-0631e63cca5c", "X6oZ1wBCflGS6ncuBufzG3lctaO8P03U"),
    ("2da5995d-9a11-442d-9160-ca2cd8dba69d", "xWJ5LHgWW3anc35kQjnS0WpRDGVVw3hR")
]
#
# The script automatically rotates through them if one fails.



if (
    not CREDENTIALS
    or CREDENTIALS[0][0].startswith("YOUR_CLIENT_ID")
):

    raise RuntimeError(
        "Sentinel Hub credentials have not been configured."
    )


cred_idx = 0


def get_config(
    idx
):

    config = SHConfig()

    config.sh_client_id = (
        CREDENTIALS[idx][0]
    )

    config.sh_client_secret = (
        CREDENTIALS[idx][1]
    )

    return config


# ================================================================
# 4. AOI
# ================================================================

BOUNDARY_URL = (
    "https://raw.githubusercontent.com/"
    "johan/world.geo.json/master/countries/"
    "SLV.geo.json"
)


print()
print("=" * 80)
print(
    "EL SALVADOR DROUGHT MONITORING"
)
print("=" * 80)


print()
print(
    "Downloading El Salvador national boundary..."
)


response = requests.get(
    BOUNDARY_URL,
    timeout=120
)

response.raise_for_status()


boundary_geojson = (
    response.json()
)


aoi_gdf = (
    gpd.GeoDataFrame.from_features(
        boundary_geojson["features"],
        crs="EPSG:4326"
    )
)


# Dissolve everything into one national polygon.

aoi_gdf = gpd.GeoDataFrame(
    geometry=[
        aoi_gdf.geometry.union_all()
    ],
    crs="EPSG:4326"
)


lon_min, lat_min, lon_max, lat_max = (
    aoi_gdf.total_bounds
)


print()
print(
    "El Salvador bounds:"
)

print(
    f"  {lon_min:.6f}, "
    f"{lat_min:.6f}, "
    f"{lon_max:.6f}, "
    f"{lat_max:.6f}"
)


aoi = BBox(
    bbox=(
        lon_min,
        lat_min,
        lon_max,
        lat_max
    ),
    crs=CRS.WGS84
)


# ================================================================
# 5. QUARTERS
# ================================================================

QUARTERS = {

    "Q1_May_Jul": [
        5,
        6,
        7
    ],

    "Q2_Aug_Oct": [
        8,
        9,
        10
    ],

    "Q3_Nov_Jan": [
        11,
        12,
        1
    ],

    "Q4_Feb_Apr": [
        2,
        3,
        4
    ]

}


# ================================================================
# 6. DETERMINE LATEST COMPLETE QUARTER
# ================================================================

TODAY = datetime.now()

current_month = TODAY.month
current_year = TODAY.year


if current_month in [5, 6, 7]:

    target_quarter_name = (
        "Q4_Feb_Apr"
    )

    target_hydro_year = (
        current_year - 1
    )


elif current_month in [8, 9, 10]:

    target_quarter_name = (
        "Q1_May_Jul"
    )

    target_hydro_year = (
        current_year
    )


elif current_month in [11, 12, 1]:

    target_quarter_name = (
        "Q2_Aug_Oct"
    )

    target_hydro_year = (
        current_year
    )


elif current_month in [2, 3, 4]:

    target_quarter_name = (
        "Q3_Nov_Jan"
    )

    target_hydro_year = (
        current_year - 1
    )


else:

    raise RuntimeError(
        "Invalid month."
    )


YEARS = [

    target_hydro_year - 2,

    target_hydro_year - 1,

    target_hydro_year

]


TARGET_QUARTER_NAME = (
    target_quarter_name
)


print()
print(
    "Running date:",
    TODAY.strftime("%Y-%m-%d")
)

print(
    "Hydrological years:",
    YEARS
)

print(
    "Target quarter:",
    TARGET_QUARTER_NAME
)

print(
    "Months:",
    QUARTERS[
        TARGET_QUARTER_NAME
    ]
)


# ================================================================
# 7. SENTINEL-2 PROCESSING PARAMETERS
# ================================================================

RESOLUTION = 90

tile_px = 500

tile_size_deg = (
    tile_px
    * RESOLUTION
    / 111320
)


# ================================================================
# 8. SPLIT EL SALVADOR INTO TILES
# ================================================================

splitter = BBoxSplitter(

    [aoi],

    CRS.WGS84,

    split_shape=(

        math.ceil(
            (lat_max - lat_min)
            /
            tile_size_deg
        ),

        math.ceil(
            (lon_max - lon_min)
            /
            tile_size_deg
        )

    )

)


tiles = (
    splitter.get_bbox_list()
)


print()
print(
    f"Generated {len(tiles)} Sentinel Hub tiles"
)


# ================================================================
# 9. INDEX DEFINITIONS
# ================================================================

INDEX_NAMES = [

    "NDDI",

    "MNDWI",

    "NMDI",

    "NDII",

    "NDVI"

]


# ================================================================
# 10. SENTINEL-2 EVALSCRIPT
# ================================================================

evalscript = """
//VERSION=3

function setup() {

    return {

        input: [

            "B03",
            "B04",
            "B08",
            "B11",
            "B12",
            "SCL",
            "dataMask"

        ],

        output: {

            id: "bands",

            bands: 5,

            sampleType: "FLOAT32"

        }

    };

}


function evaluatePixel(sample) {


    // ------------------------------------------------------------
    // CLOUD / INVALID PIXEL MASK
    // ------------------------------------------------------------

    let validSCL =

        sample.SCL === 4 ||

        sample.SCL === 5 ||

        sample.SCL === 6 ||

        sample.SCL === 7;


    let validPixel =

        sample.dataMask === 1 &&

        validSCL;


    if (!validPixel) {

        return [

            NaN,
            NaN,
            NaN,
            NaN,
            NaN

        ];

    }


    // ------------------------------------------------------------
    // SPECTRAL BANDS
    // ------------------------------------------------------------

    let GREEN = sample.B03;

    let RED = sample.B04;

    let NIR = sample.B08;

    let SWIR1 = sample.B11;

    let SWIR2 = sample.B12;


    let eps = 0.000001;


    // ------------------------------------------------------------
    // NDVI
    // ------------------------------------------------------------

    let NDVI =

        (NIR - RED) /

        (NIR + RED + eps);


    // ------------------------------------------------------------
    // NDII
    // ------------------------------------------------------------

    let NDII =

        (NIR - SWIR1) /

        (NIR + SWIR1 + eps);


    // ------------------------------------------------------------
    // NDWI USED BY NDDI
    // ------------------------------------------------------------

    let NDWI =

        (NIR - SWIR1) /

        (NIR + SWIR1 + eps);


    // ------------------------------------------------------------
    // NDDI
    // ------------------------------------------------------------

    let NDDI_denominator =

        NDVI + NDWI;


    let NDDI = NaN;


    if (

        Math.abs(
            NDDI_denominator
        ) > 0.01

    ) {

        NDDI =

            (NDVI - NDWI) /

            NDDI_denominator;

    }


    // ------------------------------------------------------------
    // MNDWI
    // ------------------------------------------------------------

    let MNDWI =

        (GREEN - SWIR1) /

        (GREEN + SWIR1 + eps);


    // ------------------------------------------------------------
    // NMDI
    // ------------------------------------------------------------

    let SWIR_difference =

        SWIR1 - SWIR2;


    let NMDI =

        (NIR - SWIR_difference) /

        (
            NIR
            +
            SWIR_difference
            +
            eps
        );


    // ------------------------------------------------------------
    // OUTPUT
    // ------------------------------------------------------------

    return [

        NDDI,

        MNDWI,

        NMDI,

        NDII,

        NDVI

    ];

}
"""


# ================================================================
# 11. FILE FUNCTIONS
# ================================================================

def composite_path(
    year,
    quarter_name
):

    return os.path.join(

        QUARTERLY_DIR,

        f"composite_{year}_"
        f"{quarter_name}.tif"

    )


def composite_exists(
    year,
    quarter_name
):

    path = composite_path(
        year,
        quarter_name
    )

    return (

        os.path.isfile(path)

        and

        os.path.getsize(path) > 0

    )


def anomaly_path(
    year,
    quarter_name
):

    return os.path.join(

        ANOMALY_DIR,

        f"anomaly_{year}_"
        f"{quarter_name}.tif"

    )


def get_month_dates(
    year,
    month
):

    last_day = (
        calendar.monthrange(
            year,
            month
        )[1]
    )


    start = datetime(
        year,
        month,
        1
    )


    end = datetime(
        year,
        month,
        last_day
    )


    return start, end


# ================================================================
# 12. SENTINEL HUB TILE REQUEST
# ================================================================

def request_tile(tile, start, end):

    global cred_idx

    attempted = 0

    while attempted < len(CREDENTIALS):

        current_cred = cred_idx

        try:

            print(
                f"\nUsing Sentinel Hub credential "
                f"{current_cred + 1}/{len(CREDENTIALS)}"
            )

            config = get_config(current_cred)

            request = SentinelHubRequest(
                evalscript=evalscript,
                input_data=[
                    SentinelHubRequest.input_data(
                        data_collection=DataCollection.SENTINEL2_L2A,
                        time_interval=(start, end)
                    )
                ],
                responses=[
                    SentinelHubRequest.output_response(
                        "bands",
                        MimeType.TIFF
                    )
                ],
                bbox=tile,
                size=(tile_px, tile_px),
                config=config
            )

            data = request.get_data()[0]

            return data

        except Exception as e:

            print()
            print(
                f"Credential {current_cred + 1} failed."
            )
            print(
                str(e)[:500]
            )

            cred_idx = (
                cred_idx + 1
            ) % len(CREDENTIALS)

            attempted += 1

    raise RuntimeError(
        "All Sentinel Hub credentials failed."
    )

# ================================================================
# 13. DOWNLOAD MONTHLY DATA / CREATE QUARTERLY TILES
# ================================================================

print()
print("=" * 80)
print(
    "SENTINEL-2 DOWNLOAD / TILE COMPOSITES"
)
print("=" * 80)


for hydro_year in YEARS:


    # ------------------------------------------------------------
    # IMPORTANT CACHE CHECK
    # ------------------------------------------------------------

    if composite_exists(

        hydro_year,

        TARGET_QUARTER_NAME

    ):

        print()

        composite_filename = os.path.basename(
            composite_path(
             	hydro_year,
                TARGET_QUARTER_NAME
            )
        )

        print(

            "Skipping Sentinel Hub download."

        )

        continue


    quarter_name = (
        TARGET_QUARTER_NAME
    )

    months = (
        QUARTERS[
            quarter_name
        ]
    )


    print()
    print("=" * 70)

    print(

        f"HYDROLOGICAL YEAR: "
        f"{hydro_year}"

    )

    print(

        f"QUARTER: "
        f"{quarter_name}"

    )

    print("=" * 70)


    for t_idx, tile in enumerate(

        tqdm(

            tiles,

            desc=
            f"{hydro_year} "
            f"{quarter_name}"

        )

    ):


        monthly_stack = []


        for month in months:


            # ----------------------------------------------------
            # January in Q3 belongs to following calendar year
            # ----------------------------------------------------

            if (

                quarter_name ==
                "Q3_Nov_Jan"

                and

                month == 1

            ):

                calendar_year = (
                    hydro_year + 1
                )

            else:

                calendar_year = (
                    hydro_year
                )


            start, end = (
                get_month_dates(
                    calendar_year,
                    month
                )
            )


            try:


                data = request_tile(

                    tile,

                    start.strftime(
                        "%Y-%m-%d"
                    ),

                    end.strftime(
                        "%Y-%m-%d"
                    )

                )


                data = (
                    data.astype(
                        np.float32
                    )
                )


                data[
                    ~np.isfinite(data)
                ] = np.nan


                monthly_stack.append(
                    data
                )


            except Exception as e:


                print()

                print(

                    f"Month failed: "
                    f"{calendar_year}-"
                    f"{month:02d}"

                )

                print(
                    str(e)[:500]
                )


        if not monthly_stack:

            continue


        with np.errstate(
            invalid="ignore"
        ):

            composite = (
                np.nanmean(

                    np.stack(

                        monthly_stack,

                        axis=0

                    ),

                    axis=0

                )
            )


        tile_path = os.path.join(

            QUARTERLY_TILE_DIR,

            f"{quarter_name}_"
            f"{hydro_year}_"
            f"tile{t_idx}.tif"

        )


        with rasterio.open(

            tile_path,

            "w",

            driver="GTiff",

            height=composite.shape[0],

            width=composite.shape[1],

            count=5,

            dtype="float32",

            crs="EPSG:4326",

            transform=Affine(
                (tile.max_x - tile.min_x) / composite.shape[1],
                0.0,
                tile.min_x,
                0.0,
                -(tile.max_y - tile.min_y) / composite.shape[0],
                tile.max_y
            ),

            compress="lzw",

            BIGTIFF="YES",

            nodata=np.nan

        ) as dst:

            for band in range(5):

                dst.write(

                    composite[:, :, band],

                    band + 1

                )

        del composite

        del monthly_stack

        gc.collect()


# ================================================================
# 14. MERGE QUARTERLY TILES
# ================================================================

print()
print("=" * 80)
print("MERGING QUARTERLY TILES")
print("=" * 80)


def merge_tiles_manual(
    tile_files,
    output_path
):
    """
    Merge regularly arranged EPSG:4326 tiles without using
    rasterio.merge.merge().

    This avoids an Affine/rasterio compatibility problem in the
    ArcGIS Pro Python environment.
    """

    if not tile_files:
        raise RuntimeError(
            "No tile files were provided for merging."
        )

    print(
        f"Merging {len(tile_files)} tiles"
    )

    # ------------------------------------------------------------
    # Open all tiles
    # ------------------------------------------------------------

    srcs = [
        rasterio.open(path)
        for path in tile_files
    ]

    try:

        # --------------------------------------------------------
        # Validate common properties
        # --------------------------------------------------------

        first = srcs[0]

        width = first.width
        height = first.height
        count = first.count
        crs = first.crs

        if count != 5:
            raise RuntimeError(
                f"Expected 5 bands, found {count}."
            )

        for src in srcs:

            if src.width != width:
                raise RuntimeError(
                    "Tile widths are inconsistent."
                )

            if src.height != height:
                raise RuntimeError(
                    "Tile heights are inconsistent."
                )

            if src.count != count:
                raise RuntimeError(
                    "Tile band counts are inconsistent."
                )

            if src.crs != crs:
                raise RuntimeError(
                    "Tile CRS values are inconsistent."
                )

        # --------------------------------------------------------
        # Get tile geographic bounds
        # --------------------------------------------------------

        tile_info = []

        for src in srcs:

            bounds = src.bounds

            tile_info.append(
                {
                    "src": src,
                    "left": float(bounds.left),
                    "right": float(bounds.right),
                    "bottom": float(bounds.bottom),
                    "top": float(bounds.top)
                }
            )

        # --------------------------------------------------------
        # Determine common pixel size
        # --------------------------------------------------------

        pixel_width = (
            tile_info[0]["right"]
            - tile_info[0]["left"]
        ) / width

        pixel_height = (
            tile_info[0]["top"]
            - tile_info[0]["bottom"]
        ) / height

        # --------------------------------------------------------
        # Determine total mosaic extent
        # --------------------------------------------------------

        west = min(
            item["left"]
            for item in tile_info
        )

        east = max(
            item["right"]
            for item in tile_info
        )

        south = min(
            item["bottom"]
            for item in tile_info
        )

        north = max(
            item["top"]
            for item in tile_info
        )

        # --------------------------------------------------------
        # Calculate mosaic dimensions
        # --------------------------------------------------------

        mosaic_width = int(
            round(
                (east - west)
                / pixel_width
            )
        )

        mosaic_height = int(
            round(
                (north - south)
                / pixel_height
            )
        )

        print(
            f"Mosaic size: "
            f"{mosaic_width} x "
            f"{mosaic_height}"
        )

        print(
            f"Pixel size: "
            f"{pixel_width:.10f} x "
            f"{pixel_height:.10f} degrees"
        )

        # --------------------------------------------------------
        # Allocate mosaic
        # --------------------------------------------------------

        mosaic = np.full(
            (
                count,
                mosaic_height,
                mosaic_width
            ),
            np.nan,
            dtype=np.float32
        )

        # --------------------------------------------------------
        # Insert each tile
        # --------------------------------------------------------

        for idx, item in enumerate(
            tile_info,
            start=1
        ):

            src = item["src"]

            print(
                f"  Adding tile "
                f"{idx}/{len(tile_info)}"
            )

            data = src.read(
                out_dtype=np.float32
            )

            # ----------------------------------------------------
            # Calculate destination position
            # ----------------------------------------------------

            col_offset = int(
                round(
                    (
                        item["left"]
                        - west
                    )
                    / pixel_width
                )
            )

            row_offset = int(
                round(
                    (
                        north
                        - item["top"]
                    )
                    / pixel_height
                )
            )

            row_end = (
                row_offset
                + height
            )

            col_end = (
                col_offset
                + width
            )

            # ----------------------------------------------------
            # Validate destination
            # ----------------------------------------------------

            if row_offset < 0:
                raise RuntimeError(
                    "Calculated negative row offset."
                )

            if col_offset < 0:
                raise RuntimeError(
                    "Calculated negative column offset."
                )

            if row_end > mosaic_height:
                raise RuntimeError(
                    "Tile extends beyond mosaic "
                    "height."
                )

            if col_end > mosaic_width:
                raise RuntimeError(
                    "Tile extends beyond mosaic "
                    "width."
                )

            # ----------------------------------------------------
            # Write tile
            # ----------------------------------------------------

            mosaic[
                :,
                row_offset:row_end,
                col_offset:col_end
            ] = data

            del data

            gc.collect()

        # --------------------------------------------------------
        # Create transform manually
        # --------------------------------------------------------

        transform = Affine(
            pixel_width,
            0.0,
            west,
            0.0,
            -pixel_height,
            north
        )

        # --------------------------------------------------------
        # Write temporary mosaic
        # --------------------------------------------------------

        temp_path = (
            output_path
            + ".tmp.tif"
        )

        profile = first.profile.copy()

        profile.update(
            driver="GTiff",
            height=mosaic_height,
            width=mosaic_width,
            count=count,
            dtype="float32",
            crs=crs,
            transform=transform,
            nodata=np.nan,
            compress="lzw",
            tiled=True,
            blockxsize=256,
            blockysize=256,
            BIGTIFF="YES"
        )

        print()
        print(
            "Writing merged mosaic..."
        )

        with rasterio.open(
            temp_path,
            "w",
            **profile
        ) as dst:

            for band in range(count):

                dst.write(
                    mosaic[band],
                    band + 1
                )

        del mosaic

        gc.collect()

        # --------------------------------------------------------
        # Return temporary mosaic
        # --------------------------------------------------------

        return temp_path

    finally:

        for src in srcs:

            src.close()

        gc.collect()


# ================================================================
# 15. CROP TO ACTUAL EL SALVADOR POLYGON
# ================================================================

def crop_raster_to_aoi(

    input_path,

    output_path,

    aoi_gdf

):


    print()

    print(

        "Cropping raster to "
        "El Salvador polygon..."

    )


    with rasterio.open(
        input_path
    ) as src:


        aoi_projected = (
            aoi_gdf.to_crs(
                src.crs
            )
        )


        geometries = [

            geometry

            for geometry
            in aoi_projected.geometry

        ]


        out_image, out_transform = mask(

            src,

            geometries,

            crop=True,

            nodata=np.nan

        )


        out_meta = (
            src.meta.copy()
        )


        out_meta.update(

            driver="GTiff",

            height=
            out_image.shape[1],

            width=
            out_image.shape[2],

            transform=
            out_transform,

            nodata=np.nan,

            dtype="float32"

        )


        with rasterio.open(

            output_path,

            "w",

            **out_meta

        ) as dst:

            dst.write(

                out_image.astype(
                    np.float32
                )

            )


    print(
        "✓ Raster cropped to "
        "national polygon"
    )


def cropped_composite_path(

    year,

    quarter_name

):

    return os.path.join(

        CROPPED_QUARTERLY_DIR,

        f"composite_{year}_"
        f"{quarter_name}_aoi.tif"

    )


print()
print("=" * 80)
print(
    "CROPPING COMPOSITES TO EL SALVADOR"
)
print("=" * 80)


for hydro_year in YEARS:


    source_path = (
        composite_path(
            hydro_year,
            TARGET_QUARTER_NAME
        )
    )


    if not os.path.exists(
        source_path
    ):

        raise RuntimeError(

            f"Missing composite: "
            f"{source_path}"

        )


    output_path = (
        cropped_composite_path(

            hydro_year,

            TARGET_QUARTER_NAME

        )
    )


    if (

        os.path.exists(output_path)

        and

        os.path.getsize(
            output_path
        ) > 0

    ):

        print()

        print(

            f"✓ Existing AOI composite: "

            f"{os.path.basename(output_path)}"

        )

        continue


    crop_raster_to_aoi(

        source_path,

        output_path,

        aoi_gdf

    )


# ================================================================
# 16. FIND AVAILABLE COMPOSITES
# ================================================================

available_years = []

raster_paths = []


for hydro_year in YEARS:


    path = (
        cropped_composite_path(

            hydro_year,

            TARGET_QUARTER_NAME

        )
    )


    if os.path.exists(path):

        available_years.append(
            hydro_year
        )

        raster_paths.append(
            path
        )

        print(
            f"✓ Found: "
            f"{os.path.basename(path)}"
        )

    else:

        print(
            f"✗ Missing: "
            f"{os.path.basename(path)}"
        )


if len(raster_paths) < 2:

    raise RuntimeError(

        "At least two yearly composites "
        "are required for anomaly calculation."

    )


# ================================================================
# 17. LEAVE-ONE-YEAR-OUT ANOMALIES
# ================================================================

print()
print("=" * 80)
print(
    "CALCULATING LEAVE-ONE-YEAR-OUT ANOMALIES"
)
print("=" * 80)


datasets = [

    rasterio.open(path)

    for path in raster_paths

]


reference = datasets[0]


for ds in datasets[1:]:


    if (

        ds.width != reference.width

        or

        ds.height != reference.height

        or

        ds.transform != reference.transform

        or

        ds.crs != reference.crs

        or

        ds.count != reference.count

    ):


        for d in datasets:
            d.close()


        raise RuntimeError(

            "Quarterly composites are "
            "not spatially compatible."

        )


anomaly_datasets = []


for hydro_year in available_years:


    output_path = (
        anomaly_path(

            hydro_year,

            TARGET_QUARTER_NAME

        )
    )


    profile = (
        reference.profile.copy()
    )


    profile.update(

        driver="GTiff",

        dtype="float32",

        count=len(INDEX_NAMES),

        compress="lzw",

        tiled=True,

        blockxsize=256,

        blockysize=256,

        BIGTIFF="YES",

        nodata=np.nan

    )


    anomaly_datasets.append(

        rasterio.open(

            output_path,

            "w",

            **profile

        )

    )


for block_number, (_, window) in enumerate(

    tqdm(

        reference.block_windows(1),

        desc="Calculating anomalies"

    )

):


    for band_idx in range(

        1,

        len(INDEX_NAMES) + 1

    ):


        year_arrays = []


        for ds in datasets:


            arr = ds.read(

                band_idx,

                window=window,

                masked=True

            )


            arr = (

                arr.filled(
                    np.nan
                )
                .astype(
                    np.float32
                )

            )


            year_arrays.append(
                arr
            )


        window_stack = np.stack(

            year_arrays,

            axis=0

        )


        valid_stack = np.isfinite(

            window_stack

        )


        total_sum = np.nansum(

            window_stack,

            axis=0,

            dtype=np.float64

        )


        total_count = np.sum(

            valid_stack,

            axis=0,

            dtype=np.int16

        )


        for i, hydro_year in enumerate(

            available_years

        ):


            current = (
                window_stack[i]
            )


            current_valid = (
                valid_stack[i]
            )


            other_sum = (

                total_sum

                -

                np.where(

                    current_valid,

                    current,

                    0

                )

            )


            other_count = (

                total_count

                -

                current_valid.astype(
                    np.int16
                )

            )


            baseline = np.full(

                current.shape,

                np.nan,

                dtype=np.float32

            )


            baseline_valid = (

                other_count > 0

            )


            baseline[
                baseline_valid
            ] = (

                other_sum[
                    baseline_valid
                ]

                /

                other_count[
                    baseline_valid
                ]

            ).astype(
                np.float32
            )


            anomaly = np.full(

                current.shape,

                np.nan,

                dtype=np.float32

            )


            anomaly_valid = (

                current_valid

                &

                baseline_valid

            )


            anomaly[
                anomaly_valid
            ] = (

                current[
                    anomaly_valid
                ]

                -

                baseline[
                    anomaly_valid
                ]

            )


            anomaly_datasets[i].write(

                anomaly,

                band_idx,

                window=window

            )


            del current
            del current_valid
            del other_sum
            del other_count
            del baseline
            del baseline_valid
            del anomaly
            del anomaly_valid


        del year_arrays
        del window_stack
        del valid_stack
        del total_sum
        del total_count


    if block_number % 50 == 0:

        gc.collect()


for dst in anomaly_datasets:

    dst.close()


for ds in datasets:

    ds.close()


gc.collect()


print()
print(
    "✓ Anomaly calculation finished"
)


# ================================================================
# 18. VERIFY ANOMALIES
# ================================================================

print()
print("=" * 80)
print(
    "VERIFYING ANOMALY OUTPUTS"
)
print("=" * 80)


for hydro_year in available_years:


    path = (
        anomaly_path(

            hydro_year,

            TARGET_QUARTER_NAME

        )
    )


    print()

    print(
        os.path.basename(path)
    )


    with rasterio.open(path) as src:


        for band_idx, index_name in enumerate(

            INDEX_NAMES,

            start=1

        ):


            finite_pixels = 0

            total_pixels = 0


            for _, window in src.block_windows(

                band_idx

            ):


                data = (

                    src.read(

                        band_idx,

                        window=window,

                        masked=True

                    )

                    .filled(np.nan)

                )


                finite_pixels += (

                    np.isfinite(
                        data
                    ).sum()

                )


                total_pixels += (
                    data.size
                )


                del data


            percentage = (

                100.0

                *

                finite_pixels

                /

                total_pixels

                if total_pixels

                else 0

            )


            print(

                f"  {index_name}: "

                f"{finite_pixels:,} / "

                f"{total_pixels:,} "

                f"({percentage:.2f}%)"

            )


# ================================================================
# 19. 5-CLASS POLYGONIZATION
# ================================================================

LOWER_PERCENTILE = 2
UPPER_PERCENTILE = 98

N_CLASSES = 5


CLASS_NAMES = {

    1: "Very Low",

    2: "Low",

    3: "Normal",

    4: "High",

    5: "Very High"

}


PERCENTILE_SAMPLE_SIZE = 2_000_000

WRITE_BATCH_SIZE = 2000

MIN_POLYGON_AREA_M2 = 14400

SIMPLIFY_TOLERANCE = 0.0001

PRESERVE_TOPOLOGY = True


def manual_window_transform(

    transform,

    window

):


    a = float(
        transform.a
    )

    b = float(
        transform.b
    )

    c = float(
        transform.c
    )

    d = float(
        transform.d
    )

    e = float(
        transform.e
    )

    f = float(
        transform.f
    )


    col_off = float(
        window.col_off
    )

    row_off = float(
        window.row_off
    )


    new_c = (

        c

        +

        col_off * a

        +

        row_off * b

    )


    new_f = (

        f

        +

        col_off * d

        +

        row_off * e

    )


    return Affine(

        a,
        b,
        new_c,
        d,
        e,
        new_f

    )


def get_percentile_sample(

    raster_path,

    band_idx,

    max_sample=PERCENTILE_SAMPLE_SIZE

):


    rng = np.random.default_rng(
        42
    )


    sample = np.empty(

        0,

        dtype=np.float32

    )


    with rasterio.open(

        raster_path

    ) as src:


        for _, window in src.block_windows(

            band_idx

        ):


            data = (

                src.read(

                    band_idx,

                    window=window,

                    masked=True

                )

                .filled(np.nan)

            )


            values = (

                data[
                    np.isfinite(data)
                ]

                .astype(
                    np.float32
                )

            )


            if values.size == 0:

                continue


            remaining = (

                max_sample

                -

                sample.size

            )


            if remaining <= 0:

                candidate = values[

                    rng.integers(

                        0,

                        values.size,

                        size=min(

                            values.size,

                            10000

                        )

                    )

                ]


                replace_count = min(

                    candidate.size,

                    sample.size

                )


                if replace_count > 0:


                    positions = rng.integers(

                        0,

                        sample.size,

                        size=replace_count

                    )


                    sample[

                        positions

                    ] = candidate[

                        :replace_count

                    ]


                continue


            if values.size <= remaining:

                sample = np.concatenate(

                    [

                        sample,

                        values

                    ]

                )

            else:

                chosen = rng.choice(

                    values,

                    size=remaining,

                    replace=False

                )


                sample = np.concatenate(

                    [

                        sample,

                        chosen

                    ]

                )


            del data
            del values


    return sample


def classify_block(

    data,

    p2,

    p98

):


    result = np.zeros(

        data.shape,

        dtype=np.uint8

    )


    valid = (

        np.isfinite(data)

        &

        (data >= p2)

        &

        (data <= p98)

    )


    if not np.any(valid):

        return result


    internal_edges = np.linspace(

        p2,

        p98,

        N_CLASSES + 1,

        dtype=np.float32

    )[1:-1]


    result[valid] = (

        np.digitize(

            data[valid],

            internal_edges,

            right=False

        )

        +

        1

    ).astype(
        np.uint8
    )


    return result


latest_year = max(
    available_years
)


latest_anomaly_path = (
    anomaly_path(

        latest_year,

        TARGET_QUARTER_NAME

    )
)


print()
print("=" * 80)
print(

    f"POLYGONIZING "
    f"{latest_year} "
    f"{TARGET_QUARTER_NAME}"

)

print("=" * 80)


polygonized_outputs = []


for band_idx, index_name in enumerate(

    INDEX_NAMES,

    start=1

):


    print()
    print("-" * 70)

    print(
        f"INDEX: {index_name}"
    )

    print("-" * 70)


    sample = get_percentile_sample(

        latest_anomaly_path,

        band_idx

    )


    if sample.size == 0:

        raise RuntimeError(

            f"No finite anomaly values "
            f"for {index_name}"

        )


    p2 = float(

        np.percentile(

            sample,

            LOWER_PERCENTILE

        )

    )


    p98 = float(

        np.percentile(

            sample,

            UPPER_PERCENTILE

        )

    )


    print(
        f"P2:  {p2:.6f}"
    )

    print(
        f"P98: {p98:.6f}"
    )


    output_name = (

        f"anomaly_{latest_year}_"

        f"{TARGET_QUARTER_NAME}_"

        f"{index_name}_classified.gpkg"

    )


    output_path = os.path.join(

        POLYGON_DIR,

        output_name

    )


    if os.path.exists(
        output_path
    ):

        os.remove(
            output_path
        )


    schema = {

        "geometry": "Polygon",

        "properties": {

            "class_id": "int",

            "class_name": "str",

            "min_value": "float",

            "max_value": "float",

            "index": "str"

        }

    }


    with rasterio.open(

        latest_anomaly_path

    ) as src:


        transform = src.transform

        crs = src.crs


        with fiona.open(

            output_path,

            "w",

            driver="GPKG",

            layer="polygons",

            schema=schema,

            crs=crs.to_string()

        ) as dst:


            batch = []


            for block_number, (_, window) in enumerate(

                tqdm(

                    src.block_windows(
                        band_idx
                    ),

                    desc=
                    f"Polygonizing {index_name}"

                )

            ):


                data = (

                    src.read(

                        band_idx,

                        window=window,

                        masked=True

                    )

                    .filled(np.nan)

                )


                classified = classify_block(

                    data,

                    p2,

                    p98

                )


                block_transform = (

                    manual_window_transform(

                        transform,

                        window

                    )

                )


                for geom, value in shapes(

                    classified,

                    mask=(
                        classified > 0
                    ),

                    transform=block_transform

                ):


                    class_id = int(
                        value
                    )


                    if class_id < 1:

                        continue


                    polygon = shape(
                        geom
                    )


                    if polygon.is_empty:

                        continue


                    if not polygon.is_valid:

                        polygon = (
                            polygon.buffer(0)
                        )


                    if polygon.is_empty:

                        continue


                    # ------------------------------------------------
                    # Area filtering
                    # ------------------------------------------------

                    polygon_projected = (

                        gpd.GeoSeries(

                            [polygon],

                            crs=crs

                        )

                        .to_crs(

                            "EPSG:32616"

                        )

                    )


                    area_m2 = float(

                        polygon_projected.area.iloc[0]

                    )


                    if (

                        area_m2
                        <
                        MIN_POLYGON_AREA_M2

                    ):

                        continue


                    if SIMPLIFY_TOLERANCE > 0:

                        polygon = polygon.simplify(

                            SIMPLIFY_TOLERANCE,

                            preserve_topology=
                            PRESERVE_TOPOLOGY

                        )


                    batch.append({

                        "geometry":
                            mapping(
                                polygon
                            ),

                        "properties": {

                            "class_id":
                                class_id,

                            "class_name":
                                CLASS_NAMES[
                                    class_id
                                ],

                            "min_value":
                                p2,

                            "max_value":
                                p98,

                            "index":
                                index_name

                        }

                    })


                    if (

                        len(batch)
                        >=
                        WRITE_BATCH_SIZE

                    ):

                        for feature in batch:

                            dst.write(
                                feature
                            )

                        batch.clear()


                del data

                del classified


                if block_number % 50 == 0:

                    gc.collect()


            for feature in batch:

                dst.write(
                    feature
                )


    polygonized_outputs.append(
        output_path
    )


    del sample

    gc.collect()


    print()

    print(
        f"✓ Saved: "
        f"{output_path}"
    )


# ================================================================
# 20. WEB OPTIMIZATION
# ================================================================

METRIC_CRS = "EPSG:32616"

WEB_CRS = "EPSG:4326"

MIN_AREA_M2 = 5000.0

SIMPLIFY_TOLERANCE_M = 20.0


def clean_geometry(
    gdf
):


    gdf = gdf[
        gdf.geometry.notna()
    ].copy()


    gdf = gdf[
        ~gdf.geometry.is_empty
    ].copy()


    invalid = (
        ~gdf.geometry.is_valid
    )


    if invalid.any():

        gdf.loc[
            invalid,
            "geometry"
        ] = (

            gdf.loc[
                invalid,
                "geometry"
            ]

            .make_valid()

        )


    gdf = gdf[
        gdf.geometry.notna()
    ].copy()


    gdf = gdf[
        ~gdf.geometry.is_empty
    ].copy()


    return gdf


print()
print("=" * 80)
print(
    "CREATING WEB-OPTIMIZED GEOJSON"
)
print("=" * 80)


web_outputs = []


for gpkg_path in polygonized_outputs:


    filename = os.path.basename(
        gpkg_path
    )


    match = re.search(

        r"_([A-Z]+)_classified\.gpkg$",

        filename

    )


    if not match:

        print(
            f"Skipping: {filename}"
        )

        continue


    index_name = (
        match.group(1)
    )


    geojson_name = (

        os.path.splitext(
            filename
        )[0]

        .replace(
            "_classified",
            ""
        )

        +

        ".geojson"

    )


    geojson_path = os.path.join(

        LATEST_VECTOR_DIR,

        geojson_name

    )


    print()

    print(
        f"Optimizing {index_name}"
    )


    gdf = gpd.read_file(

        gpkg_path,

        columns=[
            "class_id"
        ]

    )


    if gdf.empty:

        print(
            "No polygons."
        )

        continue


    gdf = clean_geometry(
        gdf
    )


    gdf = gdf.to_crs(
        METRIC_CRS
    )


    gdf["area_m2"] = (
        gdf.geometry.area
    )


    gdf = gdf[
        gdf["area_m2"]
        >=
        MIN_AREA_M2
    ].copy()


    if gdf.empty:

        print(
            "No polygons remain "
            "after area filtering."
        )

        continue


    gdf = gdf.dissolve(

        by="class_id",

        as_index=False

    )


    gdf = gdf.explode(

        index_parts=False

    ).reset_index(

        drop=True

    )


    gdf["area_m2"] = (
        gdf.geometry.area
    )


    gdf = gdf[
        gdf["area_m2"]
        >=
        MIN_AREA_M2
    ].copy()


    gdf.geometry = (
        gdf.geometry.simplify(

            SIMPLIFY_TOLERANCE_M,

            preserve_topology=True

        )
    )


    gdf = clean_geometry(
        gdf
    )


    gdf["class_id"] = (
        gdf["class_id"].astype(int)
    )


    gdf["class_name"] = (
        gdf["class_id"].map(
            CLASS_NAMES
        )
    )


    gdf["index"] = (
        index_name
    )


    gdf["min_value"] = (
        np.nan
    )


    gdf["max_value"] = (
        np.nan
    )


    gdf = gdf.to_crs(
        WEB_CRS
    )


    gdf.to_file(

        geojson_path,

        driver="GeoJSON"

    )


    web_outputs.append(
        geojson_path
    )


    print(
        f"✓ GeoJSON: "
        f"{geojson_path}"
    )


    gzip_path = (
        geojson_path +
        ".gz"
    )


    with open(

        geojson_path,

        "rb"

    ) as src:

        with gzip.open(

            gzip_path,

            "wb"

        ) as dst:

            shutil.copyfileobj(

                src,

                dst

            )


    print(
        f"✓ Gzip: "
        f"{gzip_path}"
    )


# ================================================================
# 21. COPY SIMPLE LATEST VECTOR NAMES
# ================================================================

for geojson_path in web_outputs:


    basename = os.path.basename(
        geojson_path
    )


    match = re.search(

        r"_([A-Z]+)\.geojson$",

        basename

    )


    if not match:

        continue


    index_name = (
        match.group(1)
    )


    destination = os.path.join(

        LATEST_VECTOR_DIR,

        f"{index_name}.geojson"

    )


    shutil.copy2(

        geojson_path,

        destination

    )


    print(
        f"✓ Latest vector: "
        f"{destination}"
    )


# ================================================================
# 22. COPY LATEST ANOMALY RASTERS
# ================================================================

for hydro_year in available_years:


    if hydro_year != latest_year:

        continue


    source = anomaly_path(

        hydro_year,

        TARGET_QUARTER_NAME

    )


    if os.path.exists(source):


        destination = os.path.join(

            LATEST_RASTER_DIR,

            os.path.basename(
                source
            )

        )


        shutil.copy2(

            source,

            destination

        )


        print(

            f"✓ Latest anomaly raster: "
            f"{destination}"

        )


# ================================================================
# 23. METADATA
# ================================================================

metadata = {

    "generated_at":
        TODAY.isoformat(),

    "latest_complete_quarter":
        TARGET_QUARTER_NAME,

    "latest_hydrological_year":
        latest_year,

    "hydrological_years":
        available_years,

    "quarter_months":
        QUARTERS[
            TARGET_QUARTER_NAME
        ],

    "indices":
        INDEX_NAMES,

    "resolution_m":
        RESOLUTION,

    "aoi":
        "El Salvador",

    "crs":
        "EPSG:4326",

    "anomaly_method":
        "leave-one-year-out mean",

    "outlier_percentiles": [

        LOWER_PERCENTILE,

        UPPER_PERCENTILE

    ],

    "classes":
        CLASS_NAMES,

    "minimum_polygon_area_m2":
        MIN_POLYGON_AREA_M2,

    "web_minimum_polygon_area_m2":
        MIN_AREA_M2,

    "web_simplification_m":
        SIMPLIFY_TOLERANCE_M

}


metadata_path = os.path.join(

    METADATA_DIR,

    "metadata.json"

)


with open(

    metadata_path,

    "w",

    encoding="utf-8"

) as f:


    json.dump(

        metadata,

        f,

        indent=2,

        ensure_ascii=False

    )


distribution_metadata_path = os.path.join(

    DISTRIBUTION_DIR,

    "metadata.json"

)


with open(

    distribution_metadata_path,

    "w",

    encoding="utf-8"

) as f:


    json.dump(

        metadata,

        f,

        indent=2,

        ensure_ascii=False

    )


# ================================================================
# 24. FINAL SUMMARY
# ================================================================

print()
print("=" * 80)
print(
    "DROUGHT MONITORING FINISHED SUCCESSFULLY"
)
print("=" * 80)

print()

print(
    f"Running date: "
    f"{TODAY.strftime('%Y-%m-%d')}"
)

print(
    f"Latest hydrological year: "
    f"{latest_year}"
)

print(
    f"Latest quarter: "
    f"{TARGET_QUARTER_NAME}"
)

print(
    f"Years available: "
    f"{available_years}"
)

print(
    f"Indices: "
    f"{', '.join(INDEX_NAMES)}"
)

print(
    f"Resolution: "
    f"{RESOLUTION} m"
)

print(
    f"Outlier range: "
    f"P{LOWER_PERCENTILE}-P{UPPER_PERCENTILE}"
)

print(
    f"Number of classes: "
    f"{N_CLASSES}"
)

print(
    f"Polygonized products: "
    f"{len(polygonized_outputs)}"
)

print(
    f"Web products: "
    f"{len(web_outputs)}"
)

print()

print(
    "Distribution:"
)

print(
    DISTRIBUTION_DIR
)

print()

print("=" * 80)
print("DONE")
print("=" * 80)