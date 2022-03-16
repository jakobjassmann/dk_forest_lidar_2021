# Perpare and conver Cornelius's forest distrubance map to COG
# Jakob J. Assmann j.assmann@bio.au.dk 25 February 2022

# Dependencies 
library(terra)
library(sf)
library(raster)
library(rgdal)

# Load forest change_raster
disturbance_year <- rast("data/forest_change_cs/disturbance_year_denmark.tif")

# Determine change since 2015
disturbance_since_2015 <- disturbance_year > 2015

# Load projections raster as maks
gbm_biowide <- rast("data/projections/gbm_biowide/forest_quality_gbm_biowide.vrt")

# Project change raster
disturbance_since_2015 <- terra::project(disturbance_since_2015, 
                                  gbm_biowide,
                                  method = "near")

# Mask out non-forests
disturbance_since_2015 <- mask(disturbance_since_2015, gbm_biowide)

# Write to file
writeRaster(disturbance_since_2015, 
            filename = "data/forest_change_cs/disturbance_since_2015.tif")

# Generate COG (here using gdal translate and addo)
# shell(paste("C:/OSGeo4W/OSGeo4W.bat gdal_translate",
#             "F:/JakobAssmann/dk_forest_lidar_2021/data/forest_change_cs/disturbance_since_2015_epsg3857.tif",
#             "F:/JakobAssmann/dk_forest_lidar_2021/data/forest_change_cs/disturbance_since_2015_cog_epsg3857.tif",
#             "-co TILED=YES -co COMPRESS=DEFLATE "))
# shell(paste("C:/OSGeo4W/OSGeo4W.bat gdaladdo",
#             "F:/JakobAssmann/dk_forest_lidar_2021/data/forest_change_cs/disturbance_since_2015_cog_epsg3857.tif",
#             "2 4 8 16 32"))
# 
# # This did not work, hence using GDAL addo as above
gdal_utils("warp",
           source = "data/forest_change_cs/disturbance_since_2015.tif",
           destination = "data/forest_change_cs/disturbance_since_2015_cog_epsg3857.tif",
           options = c(
             "-of", "COG",
             "-co", "RESAMPLING=NEAREST",
             "-co", "TILING_SCHEME=GoogleMapsCompatible",
             "-co", "COMPRESS=DEFLATE",
             "-co", "NUM_THREADS=46"
           ))

# EOF