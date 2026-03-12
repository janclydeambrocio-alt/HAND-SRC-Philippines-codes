# HAND-SRC-Philippines-codes

> R and Python codes for HAND-SRC-Philippines reproducibility

This repository contains the Python and R scripts developed for the thesis study:

**"Deriving HAND-Based Synthetic Rating Curves to Delineate Flood Extents in Pampanga and Cagayan de Oro River Basins"**

The scripts generate **Stage-Discharge Rating Curves (SRC)** using the **Height Above Nearest Drainage (HAND)** framework, applied to flood inundation modeling in the Philippines. They implement a dual-zone hydraulic model (channel + floodplain) with varying Manning's *n* roughness coefficients derived from land cover classification.

---

## Repository Structure

```
HAND-SRC-Philippines-codes/
├── Python/
│   ├── H1_SRC_LAD_Method.py
│   ├── H2_SRC_Einstein_Horton_Method.py
│   └── H3_SRC_Sinuosity_Effect.py
├── R/
│   └── (R scripts — coming soon)
└── README.md
```

---

## Script Descriptions

### H.1 — SRC Generation using Los Angeles Method (LAD)
**File:** `H1_SRC_LAD_Method.py`

Generates Stage-Discharge Rating Curves using the **Los Angeles Drainage (LAD)** method. Manning's *n* roughness is applied as a **pixel-count-weighted average** across land cover classes for both channel and floodplain zones. Iterates through all combinations of Manning's *n* values from literature-based ranges to produce a sensitivity ensemble of discharge estimates per stage.

---

### H.2 — SRC Generation using Einstein-Horton Method
**File:** `H2_SRC_Einstein_Horton_Method.py`

Generates Stage-Discharge Rating Curves using the **Einstein-Horton composite roughness equation**:

```
n_c = [ Σ(P_i × n_i^1.5) / P_total ]^(2/3)
```

Manning's *n* is weighted by the **wetted perimeter** of each land cover class (P_i), providing a more physically grounded composite roughness estimate than simple pixel counting. Discharge is computed separately for channel and floodplain zones and summed.

---

### H.3 — SRC Generation Considering Sinuosity Effect
**File:** `H3_SRC_Sinuosity_Effect.py`

Extends the LAD method by explicitly accounting for **channel sinuosity**. A sinuosity index (SI = 1.45) is used to adjust the energy slope, and a **1.15 multiplier** is applied to the composite Manning's *n* of the channel zone to reflect increased effective roughness due to secondary flows and energy losses at bends.

```
Se_adj  = Valley_Slope / SI
n_chan  = n_weighted_avg × 1.15
```

---

## Methodology Overview

All scripts share the same core HAND-based workflow:

1. **Inundation Raster Generation** — Raster calculator computes flood depth per pixel for stage factors from 0.1m to 3.0m
2. **Pixel-to-Point Conversion** — Each inundated pixel is converted to a point with its depth value
3. **Zone Classification** — Each pixel is classified as channel or floodplain using a channel mask shapefile
4. **Land Cover Identification** — Each pixel is assigned a land cover class (Water, Rangeland, Built Area, Crops, Trees)
5. **Hydraulic Geometry** — Cross-sectional area and wetted perimeter are computed per pixel, corrected for terrain slope
6. **Manning's n Sensitivity** — All combinations of literature-based Manning's *n* ranges are tested per land cover class
7. **Discharge Computation** — Manning's equation is applied separately for channel and floodplain zones
8. **Export** — Results are saved as Excel files per stage, and optionally written into a summary workbook

---

## Requirements

### Software
- **QGIS** (with Python console or Script Runner plugin)
- Python 3.x (bundled with QGIS)

### Python Packages
```
numpy
pandas
openpyxl
itertools (built-in)
math (built-in)
```

### QGIS Modules
```
qgis.core (QgsProject, QgsVectorLayer, QgsSpatialIndex)
processing
```

---

## Input Data Required

| Input | Format | Description |
|---|---|---|
| Water Level Raster | `.tif` | Base water level raster (PA_WL layer) |
| HAND Raster | `.tif` | Height Above Nearest Drainage |
| Slope Raster | `.tif` | Terrain slope in degrees |
| Land Cover Shapefile | `.shp` | Land cover vector with `DN` field for class ID |
| Channel Mask Shapefile | `.shp` | Polygon defining the main channel boundary |
| Summary Workbook | `.xlsx` | Existing Excel file for writing results |

> **Note:** File paths and layer names shown in the scripts above are representative. The full production script includes complete path strings specific to the study area directory structure.

---

## Manning's n Ranges (by Land Cover Class)

| Land Cover | ID | Channel n Range | Floodplain n Range |
|---|---|---|---|
| Water | 1 | 0.025 – 0.080 | 0.025 – 0.080 |
| Rangeland | 11 | 0.030 – 0.070 | 0.020 – 0.040 |
| Built Area | 7 | 0.065 – 0.120 | 0.065 – 0.120 |
| Crops | 5 | 0.025 – 0.055 | 0.025 – 0.055 |
| Trees | 2 | 0.120 – 0.300 | 0.120 – 0.300 |

> Ranges are based on Reference Roughness Literature (RRL) values.

---

## How to Run

1. Open **QGIS** and load all required raster and vector layers into your project
2. Ensure layer names in QGIS match those defined in the `# Layer Names` section of each script
3. Update all file paths in the `# USER CONFIGURATION` section
4. Open the **QGIS Python Console** (`Plugins → Python Console`)
5. Click **Show Editor**, paste or open the script, and click **Run**

---

## Key Parameters

| Parameter | Value | Description |
|---|---|---|
| `pixel_res` | 5.0 m | DEM pixel resolution |

### ⚠️ Site-Specific Parameters (Must be adjusted per station)

**All** of the following parameters are **specific to Station R10.004 (Cabula)** and must be recalibrated for other stations based on field data, DEM analysis, and observed stage-discharge records:

#### Python Scripts (H1, H2, H3)

| Parameter | Cabula Value | Description |
|---|---|---|
| `reach_length` | 3984.6 m | Hydraulic reach length — derived from DEM/channel delineation |
| `valley_S` | 0.004985 | Valley slope — computed from DEM elevation difference over distance |
| `H_crit` | 1.8 m | Critical depth for floodplain activation — based on observed data |
| `STABILITY_DEPTH_LIMIT` | 0.3 m | Minimum inundation depth for floodplain inclusion |
| `SI` | 1.45 | Sinuosity Index — measured from channel centerline vs valley length |

#### R Scripts (J1, J2, J3, J4)

| Parameter | Cabula Value | Where Used | Description |
|---|---|---|---|
| `gage_min` | 2.450 m | J1, J2, J4 | Minimum observed gage height |
| `gage_max` | 3.870 m | J1, J2, J4 | Maximum observed gage height |
| `offset_val` | -1.025 m | J1, J2, J4 | Zero-flow stage offset (h0) |
| `DEPTH_MIN` | 1.5 m | J3 | Minimum depth for sensitivity analysis |
| `DEPTH_MAX` | 2.8 m | J3 | Maximum depth for sensitivity analysis |
| `coeff_A` | 17.779 | J1, J2, J3, J4 | ORC power-law coefficient |
| `power_B` | 3.43 | J1, J2, J3, J4 | ORC power-law exponent |
| `end_row` | 74594 | J1, J2, J3, J4 | Row limit of sensitivity Excel file |

> These values are derived from the **Observed Rating Curve (ORC)** fitted in **Appendix I** (`I1_Rating_Curve_Generation.R`). Run the rating curve script first to obtain the correct values for your station before running J1–J4.

---

## Output Files

- **Flood Rasters** — GeoTIFF inundation depth rasters per stage (`Inun_0.1.tif`, `Inun_0.2.tif`, ..., `Inun_3.0.tif`)
- **Excel Files** — Sensitivity results per stage (`Sensitivity_0.1.xlsx`, ..., `Sensitivity_3.0.xlsx`)
- **Summary Workbook** — Discharge values written into an existing Excel summary sheet

---

## Study Area

**Location:** Cagayan de Oro (CDO), Philippines
**Reach:** San Isidro / Cabula subcatchment
**Projection:** Philippine coordinate system

---

## Citation

These scripts were developed as part of the following thesis study:

> **Ambrocio, J.C., et al. (2026).** *Deriving HAND-Based Synthetic Rating Curves to Delineate Flood Extents in Pampanga and Cagayan de Oro River Basins.* [Unpublished thesis].

If you use these scripts in your research, please cite the repository as:

```
Ambrocio, J.C., et al. (2026). HAND-SRC-Philippines-codes: Python and R scripts
for deriving HAND-based Synthetic Rating Curves to delineate flood extents
in Pampanga and Cagayan de Oro River Basins.
GitHub. https://github.com/janclydeambrocio-alt/HAND-SRC-Philippines-codes
```

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
