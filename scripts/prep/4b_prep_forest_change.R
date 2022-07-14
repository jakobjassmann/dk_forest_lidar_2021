# Perpare and conver Cornelius's forest distrubance map to COG
# Jakob J. Assmann j.assmann@bio.au.dk 25 February 2022

# Dependencies 
library(terra)
library(sf)
library(raster)
library(rgdal)

# Load forest change_raster
disturbance_year <- rast("data/forest_change_cs/disturbance_year_1986-2020_denmark.tif")

# Determine change since 2015
disturbance_since_2015 <- disturbance_year > 2015

# Reclassify values so that 0s are NA
disturbance_since_2015 <- classify(disturbance_since_2015,
                                   matrix(c(0, NA, 
                                            1, 1),
                                          byrow = T, ncol = 2))
# Load projections raster as template
gbm_biowide <- rast("data/projections/gbm_biowide/forest_quality_gbm_biowide.vrt")

# Project change raster
disturbance_since_2015 <- terra::project(disturbance_since_2015, 
                                         gbm_biowide,
                                         method = "near")

# Load forest mask
forest_mask <- rast("data/basemap_forests/forest_mask.tif")

# Mask out non-forests
disturbance_since_2015 <- crop(disturbance_since_2015, forest_mask)
disturbance_since_2015 <- extend(disturbance_since_2015, forest_mask)
disturbance_since_2015 <- mask(disturbance_since_2015, forest_mask)

# Crop extend to match gbm_biowide
disturbance_since_2015 <- crop(disturbance_since_2015, gbm_biowide)

# Test masking procedure
test_ext <- ext(490470,
                491360, 
                6264810,
                6265670)
forest_mask_sub <- crop(forest_mask, test_ext)
disturbance_sub <- crop(disturbance_since_2015, test_ext)
plot(forest_mask_sub)
plot(disturbance_sub)
plot(mask(disturbance_sub, forest_mask_sub))

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