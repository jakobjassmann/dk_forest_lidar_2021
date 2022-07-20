# DK Forest LiDAR v0.9.1
Classifications of Denmark's forest quality using the [EcoDes-DK15 dataset](https://github.com/jakobjassmann/ecodes-dk-lidar) and other spatial data.

**Disclaimer: This project is under development and not yet peer-reviewed.**

## Project overview
- [Workflow Overview](workflow.html)

## Data description
- [Forest annotations and training data (to be updated)](training_annotations.html)
- [Predictor overview](data_overview.html)
- [Focal (window) predictor selection](focal_var_selection.html)

## Model performance
- [Gradient Boosting performance](gbm_models_performance.html)
- [Random Forest performance](ranger_models_performance.html)

## Results
- [Leaflet web app (map of projections)](data_vis.html)
- [Summary stats (area estimates)](summary_stats.html)

## Data / Outputs
- [Summary report - website snapshot v0.1.0 (2.2 MB, PDF))](Assmann_et_al-DK_Forest_Quality_Report_v0.1.0.pdf)
- [Gradient Boosting Projections v0.9.1 (41.5 MB, GeoTiff)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/forest_quality_gbm_biowide_cog_epsg3857_v0.9.1.tif)
- [Random Forest Projections v0.9.1 (40.2 MB, GeoTiff)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/forest_quality_ranger_biowide_cog_epsg3857_v0.9.1.tif)
- [Disturbance map v0.9.1 (17.8 MB, GeoTiff)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/disturbance_since_2015_cog_epsg3857_v0.1.0.tif)
- [Training Polygons v0.9.0 (115.7 MB, GeoJson)](https://dkforestlidar2022.s3.eu-central-1.amazonaws.com/training_polygons_v0.9.0.geojson)

## Auxiliary
- [Guide on how to visualise cloud optimised rasters](cog_guide.html)

---
<img src='au_logo.png' style='width: 250; height:50px;'><img src='sustainscapes_logo.png' style='width: 240; height:60px;'>

[last update: 3 March 2022]
