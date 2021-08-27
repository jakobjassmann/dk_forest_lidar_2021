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

# get list of EcoDes-DK15 tiles
ecodes_tiles <- list.files(ecodes_vars_sub[1], pattern = "tif")
ecodes_tiles <- gsub(".*([0-9]{4}_[0-9]{3}).tif", "\\1", ecodes_tiles)

# get boundaries for bornholm
bornholm_raster <- raster("data/conif_vs_broadleaf/bornholm_forest_con_vs_dec.tif")
bornholm_bounds <- as(extent(bornholm_raster), "SpatialPolygons")
crs(bornholm_bounds) <- crs(bornholm_raster)

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
  
  # Get tile boundaries
  tile_bounds <- as(extent(tile_stack), "SpatialPolygons")
  crs(tile_bounds) <- crs(tile_stack)
  
  # Load and add forest_type data (Bornholm raster if needed)
  if(is.null(intersect(tile_bounds, bornholm_bounds))){
    forest_type_cloud <- raster("data/conif_vs_broadleaf/forest_type_cloud.tif")
    forest_type_con <- raster("data/conif_vs_broadleaf/forest_type_con.tif")
    forest_type_dec <- raster("data/conif_vs_broadleaf/forest_type_dec.tif")
  } else {
    forest_type_cloud <- raster("data/conif_vs_broadleaf/bornholm_forest_type_cloud.tif")
    forest_type_con <- raster("data/conif_vs_broadleaf/bornholm_forest_type_con.tif")
    forest_type_dec <- raster("data/conif_vs_broadleaf/bornholm_forest_type_dec.tif")
    names(forest_type_cloud) <- "forest_type_cloud"
    names(forest_type_con) <- "forest_type_con"
    names(forest_type_dec) <- "forest_type_dec"
  }
  forest_type_cloud <- crop(forest_type_cloud, tile_stack)
  forest_type_con <- crop(forest_type_con, tile_stack)
  forest_type_dec <- crop(forest_type_dec, tile_stack)
  tile_stack <- stack(tile_stack, 
                      forest_type_cloud, 
                      forest_type_con, 
                      forest_type_dec)
  
  # Add plant available water
  paw_160cm <- raster("data/plant_available_water/paw_160cm.tif")
  paw_160cm_projected <- projectRaster(paw_160cm, tile_stack)
  tile_stack <- stack(tile_stack, paw_160cm_projected)
  
  # Add focal variables
  paw_160cm_focal <- raster("data/focal_variables/a_ptv_focal_3x3.tif")
  canopy_height_focal <- raster("data/focal_variables/canopy_height_focal_3x3.tif")
  normalized_zd_sd_focal <- raster("data/focal_variables/normalized_z_sd_focal_3x3.tif")
  
  paw_160cm_focal_projected <- projectRaster(paw_160cm_focal, tile_stack)
  names(paw_160cm_focal_projected) <- "paw_160cm_focal_3x3"
  
  canopy_height_focal <- crop(canopy_height_focal, tile_stack)
  names(canopy_height_focal) <- "canopy_height_focal_3x3"
  
  normalized_zd_sd_focal <- crop(normalized_zd_sd_focal, tile_stack)
  names(normalized_zd_sd_focal) <- "normalized_zd_sd_focal_3x3"
  
  tile_stack <- stack(tile_stack, paw_160cm_focal_projected, canopy_height_focal, normalized_zd_sd_focal)
  
  # Mask water from stack
  tile_stack <- mask(tile_stack, raster(paste0(inland_mask, "/inland_water_mask_", tile_id, ".tif")))
  tile_stack <- mask(tile_stack, raster(paste0(sea_mask, "/sea_mask_", tile_id, ".tif")))
  
  # Mask non forest from raster
  forest_mask <- raster("data/projections/basemap2016_forest_mask_no_agri_forests.tif")
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
# project_tile(ecodes_tiles[12], ecodes_vars_sub, inland_mask, sea_mask, gbm_fit)

# Prepare parallel environmemt
cl <- makeCluster(74)
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(gbm))
clusterExport(cl, varlist = "bornholm_bounds")
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
