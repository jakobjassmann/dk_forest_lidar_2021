# DK Forest LiDAR by pixel projections across the EcoDes dataset
# Jakob Assmann j.assmann@bio.au.dk 9 July 2021

# Dependencies
library(caret)
library(gbm)
library(raster)
library(parallel)
library(sf)

# Load model
load("data/final_gbm_model_pixel.Rda")

# Get list of ecodes variable folders
ecodes_vars <- list.dirs("data/ecodes_subset/", recursive = T)
ecodes_vars_sub <- ecodes_vars[!grepl(".*masks", ecodes_vars)][-1]
ecodes_vars_sub <- ecodes_vars_sub[!grepl(".*vegetation_proportion$", ecodes_vars_sub)]

# get mask folders
inland_mask <- ecodes_vars[grepl("inland_water", ecodes_vars)]
sea_mask <- ecodes_vars[grepl("sea", ecodes_vars)]

# get list of ecodes tiles
ecodes_tiles <- list.files(ecodes_vars_sub[1], pattern = "tif")
ecodes_tiles <- gsub(".*([0-9]{4}_[0-9]{3}).tif", "\\1", ecodes_tiles)

# Write function to carry out predictons for one tile
project_tile <- function(tile_id, 
                         ecodes_vars_sub,
                         inland_mask,
                         sea_mask,
                         gbm_fit){
  # Load rasters as stack
  tile_stack <- stack(paste0(ecodes_vars_sub, "/",
                           gsub(".*/(.*)$", "\\1", ecodes_vars_sub),
                           "_", tile_id, ".tif"))
  names(tile_stack) <- gsub("_[0-9]{4}_[0-9]{3}", "", names(tile_stack))
  
  # Mask water from raster sack
  tile_stack <- mask(tile_stack, raster(paste0(inland_mask, "/inland_water_mask_", tile_id, ".tif")))
  tile_stack <- mask(tile_stack, raster(paste0(sea_mask, "/sea_mask_", tile_id, ".tif")))
  
  # Mask non forest from raster
  forest_mask <- raster("data/projections/basemap2016_forest_mask.tif")
  forest_mask_cropped <- crop(forest_mask, tile_stack)
  tile_stack <- mask(tile_stack, forest_mask_cropped) 

  # predict raster
  predictions_raster <- predict(tile_stack, gbm_fit)
  
  # write out raster
  writeRaster(predictions_raster, 
               paste0("data/projections/by_pixel/forest_quality_by_pixel_" , tile_id, ".tif"),
               overwrite = T)
  
  # Return nothing
  return(NULL)
}

# Prepare parallel environmnet
cl <- makeCluster(54)
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(gbm))

# Run projections
parLapply(cl, ecodes_tiles,
          project_tile,
          ecodes_vars_sub = ecodes_vars_sub,
          inland_mask = inland_mask,
          sea_mask = sea_mask,
          gbm_fit = gbm_fit)

# Stop cluster
stopCluster(cl)

# Generate VRT file
gdal_utils("buildvrt",
           source = list.files("data/projections/by_pixel/",".tif$", full.names = T),
           destination = "data/projections/by_pixel/forest_quality_by_pixel.vrt")
