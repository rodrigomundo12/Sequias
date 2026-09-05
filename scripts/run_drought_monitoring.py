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
from rasterio.windows import Window
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

def crop_raster_to_aoi(input_path, output_path, aoi_gdf):

    print()
    print("Masking composite to El Salvador boundary...")

    # ------------------------------------------------------------
    # AOI geometry
    # ------------------------------------------------------------

    aoi_geom = aoi_gdf.to_crs("EPSG:4326").geometry.unary_union

    # ------------------------------------------------------------
    # Open raster
    # ------------------------------------------------------------

    with rasterio.open(input_path) as src:

        src_crs = src.crs
        width = src.width
        height = src.height
        count = src.count

        transform = src.transform

        # --------------------------------------------------------
        # Extract affine coefficients manually
        # Avoid operations that trigger the ArcGIS/affine issue
        # --------------------------------------------------------

        a = float(transform.a)
        b = float(transform.b)
        c = float(transform.c)

        d = float(transform.d)
        e = float(transform.e)
        f = float(transform.f)

        # --------------------------------------------------------
        # Calculate raster bounds manually
        # --------------------------------------------------------

        left = c
        top = f

        right = (
            c
            + a * width
            + b * height
        )

        bottom = (
            f
            + d * width
            + e * height
        )

        raster_min_x = min(left, right)
        raster_max_x = max(left, right)

        raster_min_y = min(bottom, top)
        raster_max_y = max(bottom, top)

        # --------------------------------------------------------
        # AOI bounds
        # --------------------------------------------------------

        aoi_min_x, aoi_min_y, aoi_max_x, aoi_max_y = (
            aoi_geom.bounds
        )

        print(
            f"Raster bounds: "
            f"{raster_min_x} "
            f"{raster_min_y} "
            f"{raster_max_x} "
            f"{raster_max_y}"
        )

        print(
            f"AOI bounds: "
            f"{aoi_min_x} "
            f"{aoi_min_y} "
            f"{aoi_max_x} "
            f"{aoi_max_y}"
        )

        # --------------------------------------------------------
        # Intersection
        # --------------------------------------------------------

        min_x = max(
            raster_min_x,
            aoi_min_x
        )

        max_x = min(
            raster_max_x,
            aoi_max_x
        )

        min_y = max(
            raster_min_y,
            aoi_min_y
        )

        max_y = min(
            raster_max_y,
            aoi_max_y
        )

        if (
            min_x >= max_x
            or min_y >= max_y
        ):
            raise ValueError(
                "Raster and AOI do not overlap."
            )

        # --------------------------------------------------------
        # Pixel dimensions
        # --------------------------------------------------------

        pixel_width = abs(a)
        pixel_height = abs(e)

        # --------------------------------------------------------
        # Calculate crop window
        # --------------------------------------------------------

        col_start = int(
            np.floor(
                (min_x - raster_min_x)
                / pixel_width
            )
        )

        col_end = int(
            np.ceil(
                (max_x - raster_min_x)
                / pixel_width
            )
        )

        row_start = int(
            np.floor(
                (raster_max_y - max_y)
                / pixel_height
            )
        )

        row_end = int(
            np.ceil(
                (raster_max_y - min_y)
                / pixel_height
            )
        )

        # --------------------------------------------------------
        # Clamp to raster
        # --------------------------------------------------------

        col_start = max(
            0,
            min(col_start, width)
        )

        col_end = max(
            0,
            min(col_end, width)
        )

        row_start = max(
            0,
            min(row_start, height)
        )

        row_end = max(
            0,
            min(row_end, height)
        )

        crop_width = col_end - col_start
        crop_height = row_end - row_start

        if crop_width <= 0 or crop_height <= 0:
            raise ValueError(
                "Calculated crop window is empty."
            )

        print(
            f"Crop size: "
            f"{crop_width} x {crop_height}"
        )

        # --------------------------------------------------------
        # Read crop
        # --------------------------------------------------------

        window = rasterio.windows.Window(
            col_start,
            row_start,
            crop_width,
            crop_height
        )

        data = src.read(
            window=window
        ).astype(np.float32)

        # --------------------------------------------------------
        # Calculate new transform manually
        # --------------------------------------------------------

        new_c = (
            c
            + col_start * a
            + row_start * b
        )

        new_f = (
            f
            + col_start * d
            + row_start * e
        )

        new_transform = Affine(
            a,
            b,
            new_c,
            d,
            e,
            new_f
        )

        # --------------------------------------------------------
        # Create pixel-center coordinates
        # --------------------------------------------------------

        cols = (
            np.arange(crop_width)
            + 0.5
        )

        rows = (
            np.arange(crop_height)
            + 0.5
        )

        x_coords = (
            new_c
            + cols * a
        )

        y_coords = (
            new_f
            + rows * e
        )

        xx, yy = np.meshgrid(
            x_coords,
            y_coords
        )

        # --------------------------------------------------------
        # Determine pixels inside El Salvador
        # --------------------------------------------------------

        from shapely import contains_xy

        inside = contains_xy(
            aoi_geom,
            xx,
            yy
        )

        # --------------------------------------------------------
        # Mask everything outside the polygon
        # --------------------------------------------------------

        for band in range(count):

            band_data = data[band]

            band_data[~inside] = np.nan

            data[band] = band_data

        # --------------------------------------------------------
        # Output profile
        # --------------------------------------------------------

        profile = src.profile.copy()

        profile.update(
            driver="GTiff",
            height=crop_height,
            width=crop_width,
            count=count,
            dtype="float32",
            crs=src_crs,
            transform=new_transform,
            nodata=np.nan,
            compress="lzw",
            tiled=True,
            BIGTIFF="IF_SAFER"
        )

        # --------------------------------------------------------
        # Write cropped raster
        # --------------------------------------------------------

        with rasterio.open(
            output_path,
            "w",
            **profile
        ) as dst:

            dst.write(data)

    print(
        f"OK: Cropped composite saved: "
        f"{os.path.basename(output_path)}"
    )

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

RESOLUTION = 3000

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

    This avoids the Affine/rasterio compatibility problem found
    in the ArcGIS Pro Python environment.
    """

    if not tile_files:
        raise RuntimeError(
            "No tile files were provided for merging."
        )

    print()
    print(
        f"Merging {len(tile_files)} tiles"
    )

    # ------------------------------------------------------------
    # Open all tiles
    # ------------------------------------------------------------

    srcs = []

    try:

        for path in tile_files:

            srcs.append(
                rasterio.open(path)
            )

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

            # ------------------------------------------------------------
            # Calculate raster bounds manually.
            #
            # Avoid src.bounds because the ArcGIS Pro environment has an
            # Affine compatibility issue.
            # ------------------------------------------------------------

            transform = src.transform

            a = float(transform.a)
            b = float(transform.b)
            c = float(transform.c)

            d = float(transform.d)
            e = float(transform.e)
            f = float(transform.f)

            left = c
            top = f

            right = (
                c
                + a * src.width
                + b * src.height
            )

            bottom = (
                f
                + d * src.width
                + e * src.height
            )

            bounds = (
                min(left, right),
                min(bottom, top),
                max(left, right),
                max(bottom, top)
            )

            tile_info.append(
                {
                    "src": src,
                    "left": float(bounds[0]),
                    "bottom": float(bounds[1]),
                    "right": float(bounds[2]),
                    "top": float(bounds[3]),
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
        # Check that all tiles have the same resolution
        # --------------------------------------------------------

        tolerance = 1e-10

        for item in tile_info:

            current_pixel_width = (
                item["right"]
                - item["left"]
            ) / width

            current_pixel_height = (
                item["top"]
                - item["bottom"]
            ) / height

            if (
                abs(
                    current_pixel_width
                    - pixel_width
                ) > tolerance
            ):

                raise RuntimeError(
                    "Tile pixel widths are inconsistent."
                )

            if (
                abs(
                    current_pixel_height
                    - pixel_height
                ) > tolerance
            ):

                raise RuntimeError(
                    "Tile pixel heights are inconsistent."
                )

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
            # Replace non-finite values with NaN
            # ----------------------------------------------------

            data[
                ~np.isfinite(data)
            ] = np.nan

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
                    "Tile extends beyond mosaic height."
                )

            if col_end > mosaic_width:

                raise RuntimeError(
                    "Tile extends beyond mosaic width."
                )

            # ----------------------------------------------------
            # Write tile into mosaic
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
        # Temporary mosaic
        # --------------------------------------------------------

        temp_path = (
            output_path
            + ".tmp.tif"
        )

        # --------------------------------------------------------
        # Prepare output profile
        # --------------------------------------------------------

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

        # --------------------------------------------------------
        # Write temporary mosaic
        # --------------------------------------------------------

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

        print(
            f"Temporary mosaic created:"
            f" {os.path.basename(temp_path)}"
        )

        return temp_path

    finally:

        # --------------------------------------------------------
        # Close all source files
        # --------------------------------------------------------

        for src in srcs:

            try:
                src.close()
            except Exception:
                pass

        gc.collect()


# ================================================================
# CREATE THE NATIONAL QUARTERLY COMPOSITES
# ================================================================

print()
print(
    "Searching for quarterly tiles..."
)


for hydro_year in YEARS:

    quarter_name = TARGET_QUARTER_NAME

    # ------------------------------------------------------------
    # Find tiles belonging to this hydro year and quarter
    # ------------------------------------------------------------

    tile_files = sorted(
        [
            os.path.join(
                QUARTERLY_TILE_DIR,
                filename
            )

            for filename in os.listdir(
                QUARTERLY_TILE_DIR
            )

            if (
                filename.startswith(
                    f"{quarter_name}_{hydro_year}_"
                )
                and filename.endswith(".tif")
                and "tmp" not in filename
            )
        ]
    )

    # ------------------------------------------------------------
    # Check tile availability
    # ------------------------------------------------------------

    if not tile_files:

        print()
        print(
            f"WARNING: No quarterly tiles found "
            f"for {hydro_year}."
        )

        continue

    print()
    print("=" * 70)

    print(
        f"HYDROLOGICAL YEAR: {hydro_year}"
    )

    print(
        f"QUARTER: {quarter_name}"
    )

    print(
        f"TILES FOUND: {len(tile_files)}"
    )

    print("=" * 70)

    # ------------------------------------------------------------
    # Expected national composite
    # ------------------------------------------------------------

    output_path = composite_path(
        hydro_year,
        quarter_name
    )

    # ------------------------------------------------------------
    # If composite already exists, don't merge again
    # ------------------------------------------------------------

    if os.path.exists(output_path):

        print()
        print(
            "Existing composite found:"
        )

        print(
            f"  {os.path.basename(output_path)}"
        )

        continue

    # ------------------------------------------------------------
    # Merge tiles
    # ------------------------------------------------------------

    temp_mosaic = None

    try:

        temp_mosaic = merge_tiles_manual(
            tile_files,
            output_path
        )

        # --------------------------------------------------------
        # Crop/mask to El Salvador
        # --------------------------------------------------------

        print()
        print(
            "Masking composite to El Salvador boundary..."
        )

        crop_raster_to_aoi(
            temp_mosaic,
            output_path,
            aoi_gdf
        )

        print()
        print(
            f"Composite created:"
        )

        print(
            f"  {output_path}"
        )

    finally:

        # --------------------------------------------------------
        # Remove temporary mosaic
        # --------------------------------------------------------

        if (
            temp_mosaic is not None
            and os.path.exists(temp_mosaic)
        ):

            try:

                os.remove(
                    temp_mosaic
                )

            except Exception as e:

                print(
                    f"Warning: could not remove "
                    f"temporary file: {e}"
                )

        gc.collect()


print()
print("=" * 80)
print("QUARTERLY COMPOSITE CREATION COMPLETE")
print("=" * 80)

# ================================================================
# 15. CROP QUARTERLY COMPOSITES TO EL SALVADOR
# ================================================================

print()
print("=" * 80)
print("CROPPING COMPOSITES TO EL SALVADOR")
print("=" * 80)


def crop_raster_to_aoi(input_path, output_path, aoi_gdf):

    print()
    print("Cropping raster to El Salvador polygon...")
    print("Input:", os.path.basename(input_path))

    # ------------------------------------------------------------
    # Read AOI geometry
    # ------------------------------------------------------------

    aoi_geom = aoi_gdf.to_crs("EPSG:4326").geometry.unary_union

    # ------------------------------------------------------------
    # Open raster
    # ------------------------------------------------------------

    with rasterio.open(input_path) as src:

        src_crs = src.crs

        if src_crs is None:
            raise RuntimeError(
                "Raster has no CRS: " +
                str(input_path)
            )

        # --------------------------------------------------------
        # Read raster metadata
        # --------------------------------------------------------

        width = src.width
        height = src.height
        count = src.count

        transform = src.transform

        # Extract affine parameters explicitly.
        #
        # We do NOT call:
        #
        #     ~transform
        #
        # because the installed affine package is incompatible
        # with this ArcGIS Pro environment.
        # --------------------------------------------------------

        a = float(transform.a)
        b = float(transform.b)
        c = float(transform.c)
        d = float(transform.d)
        e = float(transform.e)
        f = float(transform.f)

        # --------------------------------------------------------
        # Raster geographic extent
        # --------------------------------------------------------

        left = c
        top = f

        right = c + a * width + b * height
        bottom = f + d * width + e * height

        raster_min_x = min(left, right)
        raster_max_x = max(left, right)

        raster_min_y = min(bottom, top)
        raster_max_y = max(bottom, top)

        print(
            "Raster bounds:",
            raster_min_x,
            raster_min_y,
            raster_max_x,
            raster_max_y
        )

        # --------------------------------------------------------
        # AOI bounds
        # --------------------------------------------------------

        aoi_min_x, aoi_min_y, aoi_max_x, aoi_max_y = (
            aoi_geom.bounds
        )

        print(
            "AOI bounds:",
            aoi_min_x,
            aoi_min_y,
            aoi_max_x,
            aoi_max_y
        )

        # --------------------------------------------------------
        # Calculate intersection
        # --------------------------------------------------------

        min_x = max(
            raster_min_x,
            aoi_min_x
        )

        max_x = min(
            raster_max_x,
            aoi_max_x
        )

        min_y = max(
            raster_min_y,
            aoi_min_y
        )

        max_y = min(
            raster_max_y,
            aoi_max_y
        )

        if min_x >= max_x or min_y >= max_y:

            raise RuntimeError(
                "AOI does not intersect raster."
            )

        # --------------------------------------------------------
        # Calculate pixel dimensions
        # --------------------------------------------------------

        pixel_width = abs(a)

        pixel_height = abs(e)

        if pixel_width == 0 or pixel_height == 0:

            raise RuntimeError(
                "Invalid raster transform."
            )

        # --------------------------------------------------------
        # Convert geographic intersection to pixel window
        #
        # This is done manually to avoid rasterio's
        # geometry_window() and Affine inversion.
        # --------------------------------------------------------

        col_start = int(
            np.floor(
                (min_x - raster_min_x) /
                pixel_width
            )
        )

        col_end = int(
            np.ceil(
                (max_x - raster_min_x) /
                pixel_width
            )
        )

        row_start = int(
            np.floor(
                (raster_max_y - max_y) /
                pixel_height
            )
        )

        row_end = int(
            np.ceil(
                (raster_max_y - min_y) /
                pixel_height
            )
        )

        # --------------------------------------------------------
        # Keep window inside raster
        # --------------------------------------------------------

        col_start = max(
            0,
            min(col_start, width)
        )

        col_end = max(
            0,
            min(col_end, width)
        )

        row_start = max(
            0,
            min(row_start, height)
        )

        row_end = max(
            0,
            min(row_end, height)
        )

        if col_start >= col_end or row_start >= row_end:

            raise RuntimeError(
                "Calculated crop window is empty."
            )

        crop_width = col_end - col_start
        crop_height = row_end - row_start

        print(
            "Crop size:",
            crop_width,
            "x",
            crop_height
        )

        # --------------------------------------------------------
        # Read crop
        # --------------------------------------------------------

        data = src.read(
            window=rasterio.windows.Window(
                col_start,
                row_start,
                crop_width,
                crop_height
            )
        ).astype(np.float32)

        # --------------------------------------------------------
        # Build output transform MANUALLY
        #
        # Again, do not use src.window_transform().
        # --------------------------------------------------------

        new_c = (
            c +
            col_start * a +
            row_start * b
        )

        new_f = (
            f +
            col_start * d +
            row_start * e
        )

        new_transform = Affine(
            a,
            b,
            new_c,
            d,
            e,
            new_f
        )

        # --------------------------------------------------------
        # Create pixel-center coordinates
        # --------------------------------------------------------

        cols = (
            np.arange(crop_width)
            + 0.5
        )

        rows = (
            np.arange(crop_height)
            + 0.5
        )

        x_coords = (
            new_c +
            cols * a
        )

        y_coords = (
            new_f +
            rows * e
        )

        # --------------------------------------------------------
        # Create coordinate mesh
        # --------------------------------------------------------

        xx, yy = np.meshgrid(
            x_coords,
            y_coords
        )

        # --------------------------------------------------------
        # Create mask WITHOUT rasterio.mask.mask()
        #
        # Shapely performs the point-in-polygon test.
        # --------------------------------------------------------

        from shapely import contains_xy

        inside = contains_xy(
            aoi_geom,
            xx,
            yy
        )

        # --------------------------------------------------------
        # Apply mask
        #
        # Everything outside El Salvador becomes NaN.
        # --------------------------------------------------------

        for band in range(count):

            band_data = data[band]

            band_data[~inside] = np.nan

            data[band] = band_data

        # --------------------------------------------------------
        # Output metadata
        # --------------------------------------------------------

        profile = src.profile.copy()

        profile.update(
            driver="GTiff",
            height=crop_height,
            width=crop_width,
            count=count,
            dtype="float32",
            crs=src_crs,
            transform=new_transform,
            nodata=np.nan,
            compress="lzw",
            tiled=True,
            BIGTIFF="IF_SAFER"
        )

        # --------------------------------------------------------
        # Write output
        # --------------------------------------------------------

        with rasterio.open(
            output_path,
            "w",
            **profile
        ) as dst:

            dst.write(data)

    print(
        "OK: Cropped composite saved:",
        os.path.basename(output_path)
    )


# ================================================================
# PROCESS ALL AVAILABLE QUARTERLY COMPOSITES
# ================================================================

for hydro_year in YEARS:

    input_path = composite_path(
        hydro_year,
        TARGET_QUARTER_NAME
    )

    output_path = os.path.join(
        QUARTERLY_DIR,
        f"composite_{hydro_year}_{TARGET_QUARTER_NAME}_cropped.tif"
    )

    if not os.path.exists(input_path):

        print()
        print(
            "WARNING: Missing composite:",
            input_path
        )

        continue

    # ------------------------------------------------------------
    # Skip if already cropped
    # ------------------------------------------------------------

    if os.path.exists(output_path):

        print()
        print(
            "OK: Cropped composite already exists:",
            os.path.basename(output_path)
        )

        continue

    crop_raster_to_aoi(
        input_path,
        output_path,
        aoi_gdf
    )


print()
print("=" * 80)
print("COMPOSITE CROPPING FINISHED")
print("=" * 80)

# ================================================================
# HELPER: GET CROPPED COMPOSITE PATH
# ================================================================

def cropped_composite_path(hydro_year, quarter_name):
    """
    Return the path of the composite cropped to the
    El Salvador national polygon.
    """

    return os.path.join(
        QUARTERLY_DIR,
        f"composite_{hydro_year}_{quarter_name}_cropped.tif"
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
            f"Found: "
            f"{os.path.basename(path)}"
        )

    else:

        print(
            f"Missing: "
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
    "Anomaly calculation finished"
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
# 19. 5-CLASS RASTER CLASSIFICATION + CLEANUP + POLYGONIZATION
# ================================================================

LOWER_PERCENTILE = 1
UPPER_PERCENTILE = 99

N_CLASSES = 5

# ================================================================
# GENERAL CLASS NAMES
# ================================================================

CLASS_NAMES = {
    1: "Muy bajo",
    2: "Bajo",
    3: "Normal",
    4: "Alto",
    5: "Muy Alto"
}

# ================================================================
# NDDI CLASS NAMES
# ================================================================

NDDI_CLASS_NAMES = {
    1: "Humedo",
    2: "Humedo / vegetacion saludable",
    3: "Normal",
    4: "Sequia moderada",
    5: "Sequia severa"
}

# ================================================================
# NDDI FIXED THRESHOLDS
#
# Class 1: < -0.1
# Class 2: -0.1 to < 0
# Class 3: 0 to < 0.1
# Class 4: 0.1 to 0.3
# Class 5: > 0.3
# ================================================================

NDDI_EDGES = [
    -np.inf,
    -0.1,
    0.0,
    0.1,
    0.3,
    np.inf
]

# ================================================================
# GENERAL SETTINGS
# ================================================================

PERCENTILE_SAMPLE_SIZE = 2_000_000

MIN_POLYGON_AREA_M2 = 144000

SIMPLIFY_TOLERANCE_M = 250

METRIC_CRS = "EPSG:32616"

TEMP_CLASSIFIED_DIR = os.path.join(
    POLYGON_DIR,
    "temporary_classified"
)

os.makedirs(
    TEMP_CLASSIFIED_DIR,
    exist_ok=True
)

MIN_REGION_PIXELS = 2

POLYGONIZE_WINDOW_SIZE = 1024


# ================================================================
# GET PERCENTILE SAMPLE
# ================================================================

def get_percentile_sample(
    raster_path,
    band_idx,
    max_samples=PERCENTILE_SAMPLE_SIZE
):

    samples = []

    with rasterio.open(raster_path) as src:

        for _, window in src.block_windows(band_idx):

            data = src.read(
                band_idx,
                window=window
            ).astype(np.float32)

            valid = data[
                np.isfinite(data)
            ]

            if valid.size == 0:
                continue

            samples.append(valid)

    if not samples:

        return np.array(
            [],
            dtype=np.float32
        )

    values = np.concatenate(samples)

    # ------------------------------------------------------------
    # Randomly reduce sample if necessary
    # ------------------------------------------------------------

    if values.size > max_samples:

        rng = np.random.default_rng(42)

        indices = rng.choice(
            values.size,
            size=max_samples,
            replace=False
        )

        values = values[indices]

    return values


# ================================================================
# CLASSIFY ONE BLOCK
# ================================================================

def classify_block(
    data,
    p2,
    p98,
    index_name
):

    result = np.zeros(
        data.shape,
        dtype=np.uint8
    )

    valid = np.isfinite(data)

    if not np.any(valid):

        return result

    # ============================================================
    # NDDI
    #
    # USE FIXED THRESHOLDS
    # ============================================================

    if index_name == "NDDI":

        # --------------------------------------------------------
        # Class 1
        # NDDI < -0.1
        # --------------------------------------------------------

        result[
            valid &
            (data < -0.1)
        ] = 1

        # --------------------------------------------------------
        # Class 2
        # -0.1 <= NDDI < 0
        # --------------------------------------------------------

        result[
            valid &
            (data >= -0.1) &
            (data < 0.0)
        ] = 2

        # --------------------------------------------------------
        # Class 3
        # 0 <= NDDI < 0.1
        # --------------------------------------------------------

        result[
            valid &
            (data >= 0.0) &
            (data < 0.1)
        ] = 3

        # --------------------------------------------------------
        # Class 4
        # 0.1 <= NDDI <= 0.3
        # --------------------------------------------------------

        result[
            valid &
            (data >= 0.1) &
            (data <= 0.3)
        ] = 4

        # --------------------------------------------------------
        # Class 5
        # NDDI > 0.3
        # --------------------------------------------------------

        result[
            valid &
            (data > 0.3)
        ] = 5

        return result

    # ============================================================
    # ALL OTHER INDEXES
    #
    # KEEP EXISTING PERCENTILE-BASED CLASSIFICATION
    # ============================================================

    clipped = np.clip(
        data,
        p2,
        p98
    )

    edges = np.linspace(
        p2,
        p98,
        N_CLASSES + 1
    )

    internal_edges = edges[1:-1]

    result[valid] = (
        np.digitize(
            clipped[valid],
            internal_edges,
            right=False
        )
        + 1
    ).astype(np.uint8)

    return result


# ================================================================
# CREATE CLASSIFIED RASTER
# ================================================================

def create_classified_raster(
    input_raster,
    output_raster,
    band_idx,
    p2,
    p98,
    index_name
):

    print()
    print("-----------------------------------------------")
    print("Creating classified raster")
    print(f"Index: {index_name}")
    print(f"Band: {band_idx}")
    print("-----------------------------------------------")

    with rasterio.open(input_raster) as src:

        profile = src.profile.copy()

        profile.update(
            dtype=rasterio.uint8,
            count=1,
            nodata=0,
            compress="lzw"
        )

        with rasterio.open(
            output_raster,
            "w",
            **profile
        ) as dst:

            for _, window in src.block_windows(
                band_idx
            ):

                data = src.read(
                    band_idx,
                    window=window
                ).astype(np.float32)

                classified = classify_block(
                    data,
                    p2,
                    p98,
                    index_name
                )

                dst.write(
                    classified,
                    1,
                    window=window
                )

    print(
        f"Classified raster saved:"
    )

    print(
        output_raster
    )


# ================================================================
# CLEAN CLASSIFIED RASTER
# ================================================================

def cleanup_classified_raster(
    input_raster,
    output_raster,
    index_name
):

    print()
    print("-----------------------------------------------")
    print("Cleaning classified raster")
    print(f"Index: {index_name}")
    print("-----------------------------------------------")

    with rasterio.open(input_raster) as src:

        profile = src.profile.copy()

        profile.update(
            dtype=rasterio.uint8,
            count=1,
            nodata=0,
            compress="lzw"
        )

        with rasterio.open(
            output_raster,
            "w",
            **profile
        ) as dst:

            for _, window in src.block_windows(1):

                data = src.read(
                    1,
                    window=window
                )

                # =================================================
                # CURRENT CLEANUP
                #
                # Preserve the classified raster.
                #
                # The previous aggressive cleanup operations
                # (majority filtering / region filling) are
                # intentionally disabled because they produced
                # excessive merging of NDDI polygons.
                # =================================================

                cleaned = data.copy()

                dst.write(
                    cleaned,
                    1,
                    window=window
                )

    print(
        f"Cleaned raster saved:"
    )

    print(
        output_raster
    )


# ================================================================
# POLYGONIZE CLEANED RASTER
# ================================================================

def polygonize_cleaned_raster(
    cleaned_raster,
    output_gpkg,
    index_name,
    p2=None,
    p98=None
):

    print()
    print("===============================================")
    print("POLYGONIZING")
    print(f"Index: {index_name}")
    print("===============================================")

    all_records = []

    # ============================================================
    # OPEN CLASSIFIED RASTER
    # ============================================================

    with rasterio.open(cleaned_raster) as src:

        transform = src.transform
        crs = src.crs

        height = src.height
        width = src.width

        # ========================================================
        # PROCESS IN WINDOWS
        # ========================================================

        for row_start in range(
            0,
            height,
            POLYGONIZE_WINDOW_SIZE
        ):

            row_end = min(
                row_start +
                POLYGONIZE_WINDOW_SIZE,
                height
            )

            for col_start in range(
                0,
                width,
                POLYGONIZE_WINDOW_SIZE
            ):

                col_end = min(
                    col_start +
                    POLYGONIZE_WINDOW_SIZE,
                    width
                )

                window = Window(
                    col_start,
                    row_start,
                    col_end - col_start,
                    row_end - row_start
                )

                data = src.read(
                    1,
                    window=window
                )

                # ------------------------------------------------
                # Ignore background
                # ------------------------------------------------

                mask = data > 0

                if not np.any(mask):
                    continue

                # =================================================
                # MANUAL WINDOW TRANSFORM
                #
                # Preserves compatibility with the rasterio
                # workflow already used in this project.
                # =================================================

                window_transform = Affine(
                    transform.a,
                    transform.b,
                    transform.c +
                    col_start * transform.a +
                    row_start * transform.b,

                    transform.d,
                    transform.e,
                    transform.f +
                    col_start * transform.d +
                    row_start * transform.e
                )

                # =================================================
                # POLYGONIZE
                # =================================================

                shapes = rasterio.features.shapes(
                    data,
                    mask=mask,
                    transform=window_transform
                )

                for geom, value in shapes:

                    class_id = int(value)

                    if class_id <= 0:
                        continue

                    all_records.append(
                        {
                            "geometry": shape(geom),
                            "class_id": class_id
                        }
                    )

    # ============================================================
    # NO POLYGONS
    # ============================================================

    if not all_records:

        print(
            f"No polygons generated for {index_name}"
        )

        return

    # ============================================================
    # CREATE GEODATAFRAME
    # ============================================================

    polygons = gpd.GeoDataFrame(
        all_records,
        crs=crs
    )

    print(
        f"Initial polygon pieces: "
        f"{len(polygons)}"
    )

    # ============================================================
    # PROJECT TO METRIC CRS
    # ============================================================

    polygons = polygons.to_crs(
        METRIC_CRS
    )

    # ============================================================
    # DISSOLVE BY CLASS
    # ============================================================

    dissolved = (
        polygons
        .dissolve(
            by="class_id",
            as_index=False
        )
    )

    print(
        f"After class dissolve: "
        f"{len(dissolved)}"
    )

    # ============================================================
    # EXPLODE MULTIPART GEOMETRIES
    # ============================================================

    dissolved = (
        dissolved
        .explode(
            index_parts=False
        )
        .reset_index(drop=True)
    )

    print(
        f"After explode: "
        f"{len(dissolved)}"
    )

    # ============================================================
    # AREA
    # ============================================================

    dissolved["area_m2"] = (
        dissolved.geometry.area
    )

    # ============================================================
    # MINIMUM AREA FILTER
    # ============================================================

    before_area = len(dissolved)

    dissolved = dissolved[
        dissolved["area_m2"] >=
        MIN_POLYGON_AREA_M2
    ].copy()

    dissolved = dissolved.reset_index(
        drop=True
    )

    print(
        f"Removed by minimum area: "
        f"{before_area - len(dissolved)}"
    )

    if dissolved.empty:

        print(
            f"No polygons remain after "
            f"area filtering for {index_name}"
        )

        return

    # ============================================================
    # SIMPLIFY GEOMETRY
    # ============================================================

    if SIMPLIFY_TOLERANCE_M > 0:

        dissolved["geometry"] = (
            dissolved.geometry.simplify(
                SIMPLIFY_TOLERANCE_M,
                preserve_topology=True
            )
        )

    # ============================================================
    # REMOVE EMPTY GEOMETRIES
    # ============================================================

    dissolved = dissolved[
        dissolved.geometry.notna()
    ].copy()

    dissolved = dissolved[
        ~dissolved.geometry.is_empty
    ].copy()

    dissolved = dissolved.reset_index(
        drop=True
    )

    # ============================================================
    # CLASS METADATA
    # ============================================================

    if index_name == "NDDI":

        # --------------------------------------------------------
        # NDDI fixed names
        # --------------------------------------------------------

        dissolved["class_name"] = (
            dissolved["class_id"]
            .map(NDDI_CLASS_NAMES)
        )

        # --------------------------------------------------------
        # NDDI fixed limits
        # --------------------------------------------------------

        nddi_min = {
            1: -np.inf,
            2: -0.1,
            3: 0.0,
            4: 0.1,
            5: 0.3
        }

        nddi_max = {
            1: -0.1,
            2: 0.0,
            3: 0.1,
            4: 0.3,
            5: np.inf
        }

        dissolved["min_value"] = (
            dissolved["class_id"]
            .map(nddi_min)
        )

        dissolved["max_value"] = (
            dissolved["class_id"]
            .map(nddi_max)
        )

    else:

        # --------------------------------------------------------
        # Other indexes keep percentile classification
        # --------------------------------------------------------

        dissolved["class_name"] = (
            dissolved["class_id"]
            .map(CLASS_NAMES)
        )

        if (
            p2 is not None and
            p98 is not None
        ):

            edges = np.linspace(
                p2,
                p98,
                N_CLASSES + 1
            )

            min_values = {
                i + 1: edges[i]
                for i in range(N_CLASSES)
            }

            max_values = {
                i + 1: edges[i + 1]
                for i in range(N_CLASSES)
            }

            dissolved["min_value"] = (
                dissolved["class_id"]
                .map(min_values)
            )

            dissolved["max_value"] = (
                dissolved["class_id"]
                .map(max_values)
            )

        else:

            dissolved["min_value"] = np.nan

            dissolved["max_value"] = np.nan

    # ============================================================
    # INDEX
    # ============================================================

    dissolved["index"] = index_name

    # ============================================================
    # AREA KM2
    # ============================================================

    dissolved["area_km2"] = (
        dissolved["area_m2"] /
        1_000_000
    )

    # ============================================================
    # RETURN TO WGS84
    # ============================================================

    dissolved = dissolved.to_crs(
        "EPSG:4326"
    )

    # ============================================================
    # SAVE GPKG
    # ============================================================

    os.makedirs(
        os.path.dirname(output_gpkg),
        exist_ok=True
    )

    # ------------------------------------------------------------
    # Remove existing GPKG if it exists
    #
    # Each index has its own GPKG, so this avoids layer conflicts.
    # ------------------------------------------------------------

    if os.path.exists(output_gpkg):

        try:

            os.remove(
                output_gpkg
            )

        except PermissionError:

            print(
                "WARNING: Could not remove existing "
                f"GeoPackage: {output_gpkg}"
            )

    # ============================================================
    # WRITE
    # ============================================================

    dissolved.to_file(
        output_gpkg,
        layer=index_name,
        driver="GPKG"
    )

    print()
    print(
        f"Final polygons: "
        f"{len(dissolved)}"
    )

    print(
        f"Saved: {output_gpkg}"
    )

    # ============================================================
    # CLASS SUMMARY
    # ============================================================

    print()
    print("Class summary:")
    print("-----------------------------------------------")

    summary = (
        dissolved
        .groupby(
            ["class_id", "class_name"],
            dropna=False
        )
        .agg(
            polygons=("class_id", "size"),
            area_km2=("area_km2", "sum")
        )
        .reset_index()
        .sort_values("class_id")
    )

    for _, row in summary.iterrows():

        print(
            f"Class {int(row['class_id'])}: "
            f"{row['class_name']} | "
            f"{int(row['polygons'])} polygons | "
            f"{row['area_km2']:.2f} km²"
        )

    print("-----------------------------------------------")

    return dissolved


# ================================================================
# RUN CLASSIFICATION + CLEANUP + POLYGONIZATION
# ================================================================

print()
print("================================================")
print("STARTING 5-CLASS PROCESSING")
print("================================================")

# ================================================================
# IMPORTANT:
#
# Process ALL five indexes.
#
# Previously:
#
# TEST_INDEXES = ["NDDI"]
#
# meant that only NDDI entered the processing loop.
# ================================================================

TEST_INDEXES = INDEX_NAMES


# ================================================================
# LATEST YEAR
# ================================================================

latest_year = max(
    available_years
)

print(
    f"Latest year: {latest_year}"
)

print(
    f"Target quarter: {TARGET_QUARTER_NAME}"
)


# ================================================================
# GET THE MULTIBAND ANOMALY RASTER
#
# THIS IS THE ORIGINAL PATH LOGIC.
# ================================================================

latest_anomaly_path = (

    anomaly_path(

        latest_year,
        TARGET_QUARTER_NAME

    )

)

print()
print(
    "Anomaly raster:"
)

print(
    latest_anomaly_path
)


# ================================================================
# PROCESS EACH INDEX
# ================================================================

for band_idx, index_name in enumerate(
    INDEX_NAMES,
    start=1
):

    # ------------------------------------------------------------
    # Only process requested indexes
    # ------------------------------------------------------------

    if index_name not in TEST_INDEXES:

        continue

    print()
    print()
    print("================================================")
    print(
        f"PROCESSING INDEX: "
        f"{index_name}"
    )
    print(
        f"Band: {band_idx}"
    )
    print("================================================")

    # ============================================================
    # VERIFY BAND
    # ============================================================

    with rasterio.open(
        latest_anomaly_path
    ) as src:

        if band_idx > src.count:

            print(
                f"WARNING: Band {band_idx} "
                f"does not exist."
            )

            print(
                f"Raster contains "
                f"{src.count} bands."
            )

            continue

    # ============================================================
    # CLASSIFICATION LIMITS
    # ============================================================

    if index_name == "NDDI":

        # --------------------------------------------------------
        # NDDI does NOT use percentiles.
        # --------------------------------------------------------

        p2 = None
        p98 = None

        print()
        print(
            "NDDI uses fixed thresholds:"
        )

        print(
            "  Class 1: NDDI < -0.1"
        )

        print(
            "  Class 2: -0.1 <= NDDI < 0"
        )

        print(
            "  Class 3: 0 <= NDDI < 0.1"
        )

        print(
            "  Class 4: 0.1 <= NDDI <= 0.3"
        )

        print(
            "  Class 5: NDDI > 0.3"
        )

    else:

        # --------------------------------------------------------
        # Other indexes retain percentile classification.
        # --------------------------------------------------------

        print()
        print(
            f"Calculating "
            f"{LOWER_PERCENTILE}th and "
            f"{UPPER_PERCENTILE}th percentiles..."
        )

        sample = get_percentile_sample(
            latest_anomaly_path,
            band_idx,
            PERCENTILE_SAMPLE_SIZE
        )

        if sample.size == 0:

            print(
                f"WARNING: No valid values "
                f"for {index_name}"
            )

            continue

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
            f"{LOWER_PERCENTILE}th percentile: "
            f"{p2}"
        )

        print(
            f"{UPPER_PERCENTILE}th percentile: "
            f"{p98}"
        )

    # ============================================================
    # CLASSIFIED RASTER
    # ============================================================

    classified_path = os.path.join(
        TEMP_CLASSIFIED_DIR,
        (
            f"{latest_year}_"
            f"{TARGET_QUARTER_NAME}_"
            f"{index_name}_classified.tif"
        )
    )

    # ============================================================
    # CLEANED RASTER
    # ============================================================

    cleaned_path = os.path.join(
        TEMP_CLASSIFIED_DIR,
        (
            f"{latest_year}_"
            f"{TARGET_QUARTER_NAME}_"
            f"{index_name}_cleaned.tif"
        )
    )

    # ============================================================
    # OUTPUT GPKG
    # ============================================================

    output_gpkg = os.path.join(
        POLYGON_DIR,
        (
            f"anomaly_"
            f"{latest_year}_"
            f"{TARGET_QUARTER_NAME}_"
            f"{index_name}_classified.gpkg"
        )
    )

    # ============================================================
    # STEP 1
    # CLASSIFY
    # ============================================================

    create_classified_raster(
        input_raster=latest_anomaly_path,
        output_raster=classified_path,
        band_idx=band_idx,
        p2=p2,
        p98=p98,
        index_name=index_name
    )

    # ============================================================
    # STEP 2
    # CLEANUP
    #
    # This is now executed for EVERY index.
    # ============================================================

    cleanup_classified_raster(
        input_raster=classified_path,
        output_raster=cleaned_path,
        index_name=index_name
    )

    # ============================================================
    # STEP 3
    # POLYGONIZATION
    # ============================================================

    polygonize_cleaned_raster(
        cleaned_raster=cleaned_path,
        output_gpkg=output_gpkg,
        index_name=index_name,
        p2=p2,
        p98=p98
    )

    # ============================================================
    # CLEAN MEMORY
    # ============================================================

    gc.collect()

    print()
    print(
        f"FINISHED: {index_name}"
    )


# ================================================================
# FINAL MESSAGE
# ================================================================

print()
print()
print("================================================")
print("ALL INDEXES FINISHED")
print("================================================")

print()
print("Processed indexes:")

for index_name in TEST_INDEXES:

    print(
        f"  - {index_name}"
    )

print()
print(
    "Temporary classified rasters:"
)

print(
    TEMP_CLASSIFIED_DIR
)

print()
print(
    "Polygon outputs:"
)

print(
    POLYGON_DIR
)

print()
print("================================================")

# ================================================================
# 20. WEB OPTIMIZATION
# ================================================================

METRIC_CRS = "EPSG:32616"

WEB_CRS = "EPSG:4326"

MIN_AREA_M2 = 5000.0

SIMPLIFY_TOLERANCE_M = 250


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
        f"GeoJSON: "
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
        f"Gzip: "
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
        f"Latest vector: "
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

            f"Latest anomaly raster: "
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