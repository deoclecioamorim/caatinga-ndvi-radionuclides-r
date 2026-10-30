# Integrating remote sensing and natural radionuclides to assess environmental variability in a Brazilian desertification nucleus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23475145.svg)](https://doi.org/10.5281/zenodo.23475145)

Code and data supporting the manuscript:

> Silva et al. (2026). *Integrating remote sensing and natural radionuclides to assess environmental variability in a Brazilian desertification nucleus*. Submitted to *Remote Sensing Applications: Society and Environment*.

The study compares two Caatinga (seasonally dry tropical forest) patches on Assunção Island, within the Cabrobó desertification nucleus (Pernambuco, Brazil), by combining:

- **Bimonthly Landsat NDVI time series (2000–2022)** analysed with BFAST (mean series of each patch, the difference between patches, and all individual pixels);
- **Natural radionuclides in soil** (⁴⁰K, ²²⁶Ra, ²²⁸Ra and ²²⁸Th; 39 plots × 3 subsamples), analysed with linear mixed models (variance partitioning) and machine learning classifiers (Random Forest and Gradient Boosting Machine, with permutation importance and SHAP values).

## Repository structure

```
R/                  analysis scripts (run in numerical order)
data/raw/           input data
  ndvi_area_A/        138 bimonthly NDVI composites of Area A (GeoTIFF, Google Earth Engine export)
  ndvi_area_B/        138 bimonthly NDVI composites of Area B
  Cabrobo_CHIRPS_2000_2022/  138 bimonthly CHIRPS precipitation totals (GeoTIFF)
  soil_radionuclides.xlsx         original laboratory activity concentrations (not included; see Data)
  plot_coordinates_area_A/B.xlsx  soil sampling plot coordinates (not included; see Data)
data/processed/     outputs of the scripts (regenerated when the scripts are run)
figures/            manuscript figures (Fig5A ... FigA7)
tables/             Table 1
gee/                Google Earth Engine script used to build the NDVI composites
```

## Google Earth Engine

`gee/ndvi_bimonthly_composites.js` builds the bimonthly NDVI composites (2000–2022) from Landsat 5, 7, 8 and 9 Collection 2 Level-2 surface reflectance: QA cloud/shadow and saturation masks, NDVI per sensor and UTC day, TM/ETM+ NDVI harmonised to OLI (Roy et al., 2016), and per-pixel median of the sensor-day observations in each two-month period, on a 30 m grid (EPSG:32724). Masked pixels are exported as NoData (-9999); no gap filling is applied. The GeoTIFFs in `data/raw/ndvi_area_A` and `data/raw/ndvi_area_B` are the output of this script (scenario S1, band `NDVI_harmOLI`). To run it, paste the file into the Earth Engine Code Editor and set `CONFIG.areaAssets` to the polygons of the study areas, which are available from the corresponding author upon reasonable request.

## Scripts

All scripts are run from the repository root (e.g. open the folder as an RStudio project, or `setwd()` to it). Each script header indicates the manuscript section it supports.

| Script | Manuscript | Output |
|---|---|---|
| `00_bfast_functions.R` | helper functions (used by 05) | – |
| `00_publication_theme.R` | figure standard (used by 02–08 and 10) | – |
| `01_ndvi_extraction.R` | Section 2.2.2 | `data/processed/ndvi_bimonthly_area_A/B.xlsx` |
| `02_bfast_mean_series.R` | Sections 2.2.2 and 3.2 | Fig. 5; `bfast_mean_series_results.xlsx` |
| `03_bfast_sensitivity_h.R` | Sections 2.2.2 and 3.2 | Fig. A.1; `bfast_sensitivity_h.xlsx` |
| `04_bfast_sensitivity_k.R` | Sections 2.2.2 and 3.2 | Fig. A.2; `bfast_sensitivity_k.xlsx` |
| `05_bfast_pixel_metrics.R` | Sections 2.2.2 and 3.2 | Fig. 6; `bfast_pixel_metrics.xlsx` and `.tif` |
| `06_bfast_spatial_resampling.R` | Sections 2.2.2 and 3.2 | Figs. A.3 and A.4; `bfast_spatial_resampling.xlsx` |
| `07_radionuclide_variance_partitioning.R` | Sections 2.3.3 and 3.3 | Fig. 7; `radionuclide_variance_partitioning.xlsx` |
| `08_radionuclide_machine_learning.R` | Sections 2.3.3 and 3.3 | Fig. 8 and Fig. A.5; `ml_results.xlsx` |
| `09_table1_descriptive_statistics.R` | Section 3.3 | Table 1 |
| `10_ndvi_precipitation.R` | NDVI × precipitation | Figs. A.6 and A.7; `ndvi_precipitation_results.xlsx` |
| `11_supplementary_material.R` | Appendix A | `Supplementary_material.docx` |

Dependencies: `01_ndvi_extraction.R` before 02–05; `05` before `06`; `02` before `10`; scripts 02–08 and 10 before `11`.

### Software

Analyses were run in R 4.6.0 with: terra 1.9.27, zoo 1.8.15, bfast 1.7.2, strucchange 1.5.4, lme4 2.0.1, caret 7.0.1, ranger 0.18.0, gbm 2.3.1, fastshap 0.3.0, shapviz 0.10.3, ggplot2 4.0.3, patchwork 1.3.2, readxl, writexl, dplyr, tibble, officer, flextable, png and ragg. Missing packages are installed automatically by each script.

## Data

**NDVI composites** (`data/raw/ndvi_area_A`, `data/raw/ndvi_area_B`): bimonthly median NDVI composites (January–February, …, November–December) from Landsat 5 TM, Landsat 7 ETM+ and Landsat 8 OLI, with TM/ETM+ harmonised to OLI; 30 m resolution, WGS 84 / UTM zone 24S (EPSG:32724). File names end in `_YYYY_MM_MM.tif` (year and bimonth).

**Precipitation** (`data/raw/Cabrobo_CHIRPS_2000_2022`): bimonthly precipitation totals (mm) from CHIRPS (0.05°, WGS 84) for a region enclosing both study areas, 2000–2022. File names end in `_YYYY_Mon_Mon.tif` (year and bimonth). Input of `10_ndvi_precipitation.R`.

**Soil radionuclides** (`data/processed/soil_radionuclides.xlsx`, sheet `samples`; not included in this version): activity concentrations (Bq kg⁻¹) determined by high-resolution gamma spectrometry in soil samples from the 0–20 cm layer. This file is the input of scripts 07–09; the results of these scripts (variance components, classification performance and variable importance, Table 1) are provided in `data/processed/` and `tables/`.

| Column | Description |
|---|---|
| `area` | study patch (A or B) |
| `plot` | sampling plot (Area A: 1–21; Area B: 1–18) |
| `subsample` | subsample within the plot (1–3) |
| `ra228`, `th228`, `ra226`, `k40` | activity concentrations of ²²⁸Ra, ²²⁸Th, ²²⁶Ra and ⁴⁰K (Bq kg⁻¹) |

The radionuclide data and the soil sampling plot coordinates (`data/raw/plot_coordinates_area_A/B.xlsx`), required by scripts 05 and 07–09, are not included in this version and are available from the corresponding author upon reasonable request.

## Citation

If you use this code or data, please cite the article above (DOI to be added on publication) and this repository:

> Silva, B. F. da, Bertaco, G. P., Silva, P. S. C. da, Rodrigues, R. Z. dos S., Santos, M. da S., Pereira, E. C., Amorim, D. J., Cardoso Junior, M. J. de M., Ballester, M. V. R., Heiskanen, J., Silva, A. K. de O., & Camargo, P. B. de (2026). *Code and data for: Integrating remote sensing and natural radionuclides to assess environmental variability in a Brazilian desertification nucleus*. Zenodo. https://doi.org/10.5281/zenodo.23475145

The DOI above always resolves to the latest version of this repository.

## License

Released under the MIT License (see `LICENSE`).
