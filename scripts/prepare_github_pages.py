# ================================================================
# PREPARE GITHUB PAGES
#
# Copies the latest drought-monitoring GeoJSON files to:
#     docs/vectors/
#
# Creates:
#     docs/vectors/index.html
#
# The web map:
#   - Displays five drought indices
#   - Allows ONLY ONE index to be visible at a time
#   - Uses checkbox-style controls
#   - Starts with NDDI selected
#   - Uses GeoJSON properties:
#         class = 1...5
#         label = Very Low ... Very High
# ================================================================

from pathlib import Path
import shutil
import sys


# ================================================================
# 1. PROJECT DIRECTORIES
# ================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

SOURCE_VECTORS = (
    PROJECT_ROOT
    / "distribution"
    / "latest"
    / "vectors"
)

WEB_VECTORS = (
    PROJECT_ROOT
    / "docs"
    / "vectors"
)

SOURCE_METADATA = (
    PROJECT_ROOT
    / "metadata"
)

WEB_METADATA = (
    PROJECT_ROOT
    / "docs"
    / "metadata"
)


# ================================================================
# 2. DROUGHT INDICES
# ================================================================

INDICES = [
    "NDDI",
    "MNDWI",
    "NMDI",
    "NDII",
    "NDVI",
]

# ================================================================
# QUARTER / MONTH TITLE
# ================================================================

QUARTER_MONTHS = {
    "Q1_May_Jul": "MAYO a JULIO",
    "Q2_Aug_Oct": "AGOSTO a OCTUBRE",
    "Q3_Nov_Jan": "NOVIEMBRE a ENERO",
    "Q4_Feb_Apr": "FEBRERO a ABRIL",
}

# ================================================================
# 3. CLASSIFICATION
# ================================================================

CLASS_COLORS = {
    1: "#F4B4B4",
    2: "#FAD09E",
    3: "#FFF9A6",
    4: "#B2E2E2",
    5: "#AEC6CF"
    #1: "#8B0000",   # Very Low
    #2: "#FF8C00",   # Low
    #3: "#FFFF66",   # Normal
    #4: "#90EE90",   # High
    #5: "#006400",   # Very High
}

CLASS_NAMES = {
    1: "Very Low",
    2: "Low",
    3: "Normal",
    4: "High",
    5: "Very High",
}


# ================================================================
# 4. CHECK PROJECT DIRECTORIES
# ================================================================

print()
print("=" * 60)
print("PREPARING GITHUB PAGES")
print("=" * 60)

print()
print("Project directory:")
print(PROJECT_ROOT)

print()
print("Source vectors:")
print(SOURCE_VECTORS)

print()
print("Web vectors:")
print(WEB_VECTORS)


if not PROJECT_ROOT.exists():
    print()
    print("ERROR: Project directory does not exist.")
    sys.exit(1)


if not SOURCE_VECTORS.exists():
    print()
    print("ERROR: Source vector directory does not exist:")
    print(SOURCE_VECTORS)
    sys.exit(1)


# ================================================================
# 5. CLEAN PREVIOUS WEB VECTORS
# ================================================================

print()
print("=" * 60)
print("CLEANING PREVIOUS WEB VECTORS")
print("=" * 60)

if WEB_VECTORS.exists():
    print("Removing:")
    print(WEB_VECTORS)

    shutil.rmtree(WEB_VECTORS)

WEB_VECTORS.mkdir(
    parents=True,
    exist_ok=True
)

print("Created:")
print(WEB_VECTORS)


# ================================================================
# 6. COPY GEOJSON FILES
# ================================================================

print()
print("=" * 60)
print("COPYING GEOJSON FILES")
print("=" * 60)

missing_files = []

for index in INDICES:

    source = SOURCE_VECTORS / f"{index}.geojson"
    destination = WEB_VECTORS / f"{index}.geojson"

    print()
    print(f"{index}:")

    if not source.exists():

        print("  ERROR: Missing source file:")
        print(f"  {source}")

        missing_files.append(index)

        continue

    shutil.copy2(
        source,
        destination
    )

    print("  Source:")
    print(f"    {source}")

    print("  Destination:")
    print(f"    {destination}")

    print("  OK")


if missing_files:

    print()
    print("=" * 60)
    print("ERROR: MISSING GEOJSON FILES")
    print("=" * 60)

    for index in missing_files:
        print(f" - {index}.geojson")

    print()
    print("The web page will NOT be generated.")
    sys.exit(1)

# ================================================================
# DETECT CURRENT QUARTER
# ================================================================

current_quarter = None

if SOURCE_METADATA.exists():

    metadata_files = list(
        SOURCE_METADATA.glob("*.json")
    )

    for metadata_file in metadata_files:

        try:

            metadata_text = metadata_file.read_text(
                encoding="utf-8"
            )

            for quarter_name in QUARTER_MONTHS:

                if quarter_name in metadata_text:

                    current_quarter = quarter_name
                    break

            if current_quarter:
                break

        except Exception:
            continue


if current_quarter is None:

    print()
    print(
        "WARNING: Could not automatically detect quarter."
    )

    title_months = "PERIODO NO DISPONIBLE"

else:

    title_months = (
        QUARTER_MONTHS[current_quarter]
    )

    print()
    print(
        f"Detected quarter: {current_quarter}"
    )

    print(
        f"Map title: ANOMALIAS DE {title_months}"
    )

# ================================================================
# 7. CREATE WEB MAP HTML
# ================================================================

print()
print("=" * 60)
print("CREATING WEB MAP")
print("=" * 60)


MAP_HTML = r"""<!DOCTYPE html>
<html lang="en">

<head>

    <meta charset="UTF-8">

    <meta
        name="viewport"
        content="width=device-width, initial-scale=1.0"
    >

    <title>El Salvador Drought Monitoring</title>

    <link
        rel="stylesheet"
        href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"
    >

    <style>

        html,
        body {
            margin: 0;
            padding: 0;
            width: 100%;
            height: 100%;
        }

        #map {
            width: 100%;
            height: 100vh;
        }

        /* -------------------------------------------------------
           INDEX CONTROL
           ------------------------------------------------------- */

        .index-control {
            background: white;
            padding: 10px 12px;
            border-radius: 5px;
            box-shadow: 0 1px 5px rgba(0,0,0,0.4);
            font-family: Arial, sans-serif;
            line-height: 1.8;
        }

        .index-control-title {
            font-weight: bold;
            margin-bottom: 5px;
            font-size: 14px;
        }

        .index-option {
            display: block;
            white-space: nowrap;
            cursor: pointer;
        }

        .index-option input {
            margin-right: 6px;
            cursor: pointer;
        }

        /* -------------------------------------------------------
           LEGEND
           ------------------------------------------------------- */

        .legend {
            background: white;
            padding: 10px 12px;
            border-radius: 5px;
            box-shadow: 0 1px 5px rgba(0,0,0,0.4);
            font-family: Arial, sans-serif;
            line-height: 18px;
        }

        .legend-title {
            font-weight: bold;
            margin-bottom: 5px;
        }

        /* -------------------------------------------------------
           MAP TITLE
           ------------------------------------------------------- */

        .map-title {
            background: white;
            padding: 8px 16px;
            border-radius: 5px;
            box-shadow: 0 1px 5px rgba(0,0,0,0.4);
            font-family: Arial, sans-serif;
            font-weight: bold;
            font-size: 16px;
            text-align: center;
        }

        .legend-item {
            display: flex;
            align-items: center;
        }

        .legend-color {
            width: 18px;
            height: 18px;
            margin-right: 7px;
            border: 1px solid #555;
        }

        .leaflet-control.map-title-control {
            position: fixed;
            left: 50%;
            transform: translateX(-50%);
            top: 10px;
            margin: 0;
        }

    </style>

</head>


<body>

<div id="map"></div>


<script
    src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js">
</script>


<script>

    // ============================================================
    // CONFIGURATION
    // ============================================================

    const indices = [
        "NDDI",
        "MNDWI",
        "NMDI",
        "NDII",
        "NDVI"
    ];


    // ============================================================
    // MAP
    // ============================================================

    const map = L.map("map").setView(
        [13.7942, -88.8965],
        8
    );

    // ============================================================
    // MAP TITLE
    // ============================================================


    const mapTitle =
        L.control({
            position: "topright"
        });


    mapTitle.onAdd = function() {

        const div =
            L.DomUtil.create(
                "div",
                "map-title"
            );

        div.innerHTML =
            "ANOMALÍAS DE __TITLE_MONTHS__";

        return div;
    };

    mapTitle.addTo(map);



    // ============================================================
    // BASE MAP
    // ============================================================

    const osm = L.tileLayer(
        "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
        {
            maxZoom: 19,
            attribution:
                '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>'
        }
    ).addTo(map);


    // ============================================================
    // CLASSIFICATION
    // ============================================================

    const classColors = {

        1: "#F4B4B4",
        2: "#FAD09E",
        3: "#FFF9A6",
        4: "#B2E2E2",
        5: "#AEC6CF"

    };


    const classNames = {

        1: "Very Low",
        2: "Low",
        3: "Normal",
        4: "High",
        5: "Very High"

    };


    // ============================================================
    // STORAGE FOR GEOJSON LAYERS
    // ============================================================

    const overlayLayers = {};


    // ============================================================
    // STYLE GEOJSON
    //
    // IMPORTANT:
    // The GeoJSON properties are:
    //
    //     class
    //     label
    //
    // NOT:
    //
    //     class_id
    //     class_name
    // ============================================================

    function styleFeature(feature) {

    const properties = feature.properties || {};

    // Support both GeoJSON schemas:
    //
    // NDDI:
    //   class_id
    //
    // Other indices:
    //   class
    //
    const classValue = Number(
        properties.class ?? properties.class_id
    );

    return {

            color: "#444",

            weight: 0.4,

            opacity: 0.8,

            fillColor:
                classColors[classValue] || "#999",

            fillOpacity: 0.65

        };

    }


    // ============================================================
    // POPUP
    // ============================================================

    function onEachFeature(feature, layer) {

    const properties =
        feature.properties || {};


    // Support both GeoJSON schemas:
    //
    // NDDI:
    //   class_id
    //   class_name
    //
    // Other indices:
    //   class
    //   label
    //
    const classValue = Number(
        properties.class ??
        properties.class_id
    );


    const label =
        properties.label ??
        properties.class_name ??
        classNames[classValue] ??
        "Unknown";


        const popupHTML =

            "<div style='font-family: Arial, sans-serif;'>" +

            //"<b>Index:</b> " +
            //currentIndex +

            //"<br>" +

            "<b>Class:</b> " +
            label +

            "<br>" +

            "<b>Class value:</b> " +
            (
                Number.isFinite(classValue)
                    ? classValue
                    : "N/A"
            ) +

            "</div>";


        layer.bindPopup(
            popupHTML
        );

    }


    // ============================================================
    // CURRENTLY SELECTED INDEX
    // ============================================================

    let currentIndex = "NDDI";


    // ============================================================
    // LOAD ALL GEOJSON FILES
    // ============================================================

    async function loadAllLayers() {

        const promises = indices.map(
            async function(index) {

                const response = await fetch(
                    index + ".geojson"
                );


                if (!response.ok) {

                    throw new Error(
                        "Could not load " +
                        index +
                        ".geojson"
                    );

                }


                const data =
                    await response.json();


                const layer = L.geoJSON(
                    data,
                    {
                        style: styleFeature,
                        onEachFeature: onEachFeature
                    }
                );


                overlayLayers[index] =
                    layer;


                return index;

            }
        );


        await Promise.all(
            promises
        );


        // ========================================================
        // SHOW NDDI INITIALLY
        // ========================================================

        overlayLayers["NDDI"].addTo(map);


        currentIndex = "NDDI";


        // ========================================================
        // CREATE SINGLE-SELECTION CONTROL
        // ========================================================

        createIndexControl();

    }


    // ============================================================
    // CREATE INDEX CONTROL
    //
    // The controls look like checkboxes.
    //
    // However, JavaScript forces them to behave as a
    // SINGLE-SELECTION control.
    //
    // Therefore:
    //
    // NDDI     [x]
    // MNDWI    [ ]
    // NMDI     [ ]
    // NDII     [ ]
    // NDVI     [ ]
    //
    // Selecting another index automatically hides the
    // previously selected one.
    // ============================================================

    function createIndexControl() {

        const control =
            L.control({
                position: "topright"
            });


        control.onAdd = function() {

            const container =
                L.DomUtil.create(
                    "div",
                    "index-control"
                );


            L.DomEvent.disableClickPropagation(
                container
            );


            const title =
                L.DomUtil.create(
                    "div",
                    "index-control-title",
                    container
                );


            title.innerHTML =
                "Anomalias de indices multiespectrales";


            indices.forEach(
                function(index) {

                    const label =
                        L.DomUtil.create(
                            "label",
                            "index-option",
                            container
                        );


                    const checkbox =
                        L.DomUtil.create(
                            "input",
                            "",
                            label
                        );


                    checkbox.type =
                        "checkbox";


                    checkbox.value =
                        index;


                    checkbox.checked =
                        (
                            index ===
                            currentIndex
                        );


                    const text =
                        document.createTextNode(
                            index
                        );


                    label.appendChild(
                        text
                    );


                    // ------------------------------------------------
                    // SINGLE-SELECTION BEHAVIOR
                    // ------------------------------------------------

                    checkbox.addEventListener(
                        "change",
                        function() {

                            // If user tries to uncheck the
                            // currently active index, immediately
                            // turn it back on.
                            if (!checkbox.checked) {

                                checkbox.checked = true;

                                return;

                            }


                            // Hide every other index.
                            indices.forEach(
                                function(otherIndex) {

                                    if (otherIndex !== index) {

                                        const otherLayer =
                                            overlayLayers[otherIndex];


                                        if (
                                            otherLayer &&
                                            map.hasLayer(otherLayer)
                                        ) {

                                            map.removeLayer(
                                                otherLayer
                                            );

                                        }

                                    }

                                }
                            );


                            // Show selected index.
                            const selectedLayer =
                                overlayLayers[index];


                            if (
                                selectedLayer &&
                                !map.hasLayer(selectedLayer)
                            ) {

                                selectedLayer.addTo(map);

                            }


                            // Update current index.
                            currentIndex =
                                index;


                            // Update all checkbox states.
                            const allCheckboxes =
                                container.querySelectorAll(
                                    "input[type='checkbox']"
                                );


                            allCheckboxes.forEach(
                                function(otherCheckbox) {

                                    otherCheckbox.checked =
                                        (
                                            otherCheckbox.value ===
                                            currentIndex
                                        );

                                }
                            );

                        }
                    );

                }
            );


            return container;

        };


        control.addTo(
            map
        );

    }


    // ============================================================
    // LEGEND
    // ============================================================

    const legend =
        L.control({
            position: "bottomright"
        });


    legend.onAdd = function() {

        const div =
            L.DomUtil.create(
                "div",
                "legend"
            );


        div.innerHTML =
            "<div class='legend-title'>" +
            "Clasificacion" +
            "</div>";


        const classes = [
            [1, "Muy bajo"],
            [2, "Bajo"],
            [3, "Normal"],
            [4, "Alto"],
            [5, "Muy Alto"]
        ];


        classes.forEach(
            function(item) {

                const classValue =
                    item[0];

                const classLabel =
                    item[1];


                div.innerHTML +=

                    "<div class='legend-item'>" +

                    "<span class='legend-color' " +
                    "style='background:" +
                    classColors[classValue] +
                    ";'></span>" +

                    classLabel +

                    "</div>";

            }
        );


        return div;

    };


    legend.addTo(
        map
    );


    // ============================================================
    // LOAD DATA
    // ============================================================

    loadAllLayers()
        .catch(
            function(error) {

                console.error(
                    "Error loading drought layers:",
                    error
                );


                alert(
                    "Error loading drought data. " +
                    "Please check the browser console."
                );

            }
        );

</script>

</body>
</html>
"""

MAP_HTML = MAP_HTML.replace(
    "__TITLE_MONTHS__",
    title_months
)


# ================================================================
# 8. WRITE MAP HTML
# ================================================================

MAP_FILE = (
    WEB_VECTORS
    / "index.html"
)

MAP_FILE.write_text(
    MAP_HTML,
    encoding="utf-8"
)

print()
print("Map created:")
print(MAP_FILE)


# ================================================================
# 9. COPY METADATA
# ================================================================

print()
print("=" * 60)
print("COPYING METADATA")
print("=" * 60)


if SOURCE_METADATA.exists():

    if WEB_METADATA.exists():
        shutil.rmtree(
            WEB_METADATA
        )

    WEB_METADATA.mkdir(
        parents=True,
        exist_ok=True
    )


    for item in SOURCE_METADATA.iterdir():

        destination = WEB_METADATA / item.name

        if item.is_dir():

            shutil.copytree(
                item,
                destination
            )

        else:

            shutil.copy2(
                item,
                destination
            )

        print(
            f"Copied: {item.name}"
        )

else:

    print(
        "Metadata directory not found."
    )

    print(
        "Skipping metadata copy."
    )


# ================================================================
# 10. VERIFY OUTPUTS
# ================================================================

print()
print("=" * 60)
print("VERIFYING WEB PRODUCTS")
print("=" * 60)


verification_failed = False


for index in INDICES:

    file_path = WEB_VECTORS / f"{index}.geojson"

    if file_path.exists():

        size_mb = file_path.stat().st_size / (1024 * 1024)

        print(
            f"OK: {index}.geojson "
            f"({size_mb:.2f} MB)"
        )

    else:

        print(
            f"ERROR: {index}.geojson missing"
        )

        verification_failed = True


if MAP_FILE.exists():

    print(
        "OK: index.html"
    )

else:

    print(
        "ERROR: index.html missing"
    )

    verification_failed = True


# ================================================================
# 11. VERIFY IMPORTANT JAVASCRIPT
# ================================================================

print()
print("Checking map JavaScript...")

html_text = MAP_FILE.read_text(
    encoding="utf-8"
)


required_strings = [

    'const indices = [',

    '"NDDI"',

    '"MNDWI"',

    '"NMDI"',

    '"NDII"',

    '"NDVI"',

    'properties.class',

    'properties.label',

    'createIndexControl',

    'checkbox.type',

    'checkbox.checked',

    'currentIndex = "NDDI"',

    'overlayLayers["NDDI"].addTo(map)',

]


for required in required_strings:

    if required in html_text:

        print(
            f"  OK: {required}"
        )

    else:

        print(
            f"  ERROR: Missing: {required}"
        )

        verification_failed = True


# ================================================================
# 12. FINAL STATUS
# ================================================================

print()
print("=" * 60)

if verification_failed:

    print("PREPARATION FAILED")

    print("=" * 60)

    sys.exit(1)

else:

    print("GITHUB PAGES PREPARATION COMPLETE")

    print("=" * 60)

    print()
    print("Web vectors:")
    print(WEB_VECTORS)

    print()
    print("Map:")
    print(MAP_FILE)

    print()
    print("Available indices:")

    for index in INDICES:
        print(f"  - {index}")

    print()
    print("Initial index: NDDI")

    print()
    print(
        "Only one drought index can be displayed at a time."
    )

    print()
    print(
        "GeoJSON properties recognized:"
    )

    print(
        "  class"
    )

    print(
        "  label"
    )

    print()