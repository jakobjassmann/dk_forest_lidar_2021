# High conservation value forest in Denmark (DK Forest LiDAR v1.0.0)

**Assmann et al. (in prep) Temperate forests of high conservation value are successfully identified by satellite and LiDAR data fusion.**

Classifications of high conservation value of forests in Denmark using the [EcoDes-DK15 dataset](https://github.com/jakobjassmann/ecodes-dk-lidar) and other spatial data.

**Disclaimer: This project is currently in peer-review.**

## Results
- [Leaflet web app - best model only (map of projections)](data_vis_best.html)
- [Leaflet web app - all models (map of projections)](data_vis_all.html)
- [Summary stats (area estimates)](summary_stats.html)

## Project overview
- [Workflow Overview](workflow.html)

## Data description
- [Forest annotations and training data](training_annotations.html)
- [Predictor overview](data_overview.html)
- [Focal (window) predictor selection](focal_var_selection.html)

## Model performance
- [Gradient Boosting performance](gbm_models_performance.html)
- [Random Forest performance](ranger_models_performance.html)

## Data / Outputs
- [Summary report - website snapshot v1.0.0 (3.48 MB, PDF))](Assmann_et_al-DK_Forest_Quality_Report_v1.0.0.pdf)
- [Best model: Random Forest Projections BIOWIDE v1.0.0 (40.2 MB, GeoTiff)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/forest_quality_ranger_biowide_cog_epsg3857_v1.0.0.tif)
- [Random Forest Projections Derek's Stratification v1.0.0 (40.2 MB, GeoTiff)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/forest_quality_ranger_derek_cog_epsg3857_v1.0.0.tif)
- [Gradient Boosting Projections BIOWIDE v1.0.0 (41.5 MB, GeoTiff)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/forest_quality_gbm_biowide_cog_epsg3857_v1.0.0.tif)
- [Gradient Boosting Projections Derek's Stratification v1.0.0 (41.5 MB, GeoTiff)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/forest_quality_gbm_derekcog_epsg3857_v1.0.0.tif)
- [Disturbance map v0.9.1 (17.8 MB, GeoTiff)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/disturbance_since_2015_cog_epsg3857_v0.1.0.tif)
- [Training Polygons v0.9.0 (115.7 MB, GeoJson)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/training_polygons_v0.9.0.geojson)

## Auxiliary
- [Guide on how to visualise cloud optimised rasters](cog_guide.html)

---
<img src='au_logo.png' style='width: 250; height:50px;'><img src='sustainscapes_logo.png' style='width: 240; height:60px;'>

[last update: 9 August 2022]
