# DK Forest LiDAR by pixel projections across the EcoDes dataset
# Jakob Assmann j.assmann@bio.au.dk 9 July 2021

## 1) Data and environment prep ----

# Dependencies
library(caret)
library(gbm)
library(raster)
library(terra)
library(parallel)
library(pbapply)
library(sf)
library(rgdal)

# Load model
load("data/models/final_gbm_model_pixel_biowide.Rda")

# Rename model 
model_fit <- gbm_fit
rm(gbm_fit)

# Get list predictor variables 
pred_vars<- summary(model_fit, plotit = F)$var

# Get list of ecodes vrt files available (using dir here to speed up listing)
ecodes_dir <- "F:/JakobAssmann/EcoDes-DK15_v1.1.0/" 
ecodes_vars <- shell(paste0("dir /b /s ", 
                              gsub("/", "\\\\", ecodes_dir),
                              "*.vrt"),
                       intern = T) %>%
  gsub("\\\\", "/", .) %>%
  gsub("(.*/).*vrt$", "\\1", .)

# Filter out ecodes predictor files
ecodes_preds <- sapply(ecodes_vars, 
                       function(x) gsub(".*/(.*)/", "\\1", x) %in% pred_vars) %>%
  ecodes_vars[.]

# Get list of remaining predictors
pred_files <- list.files("data/predictor_data/", "tif", recursive = T, 
                         full.names = T)
pred_files <- sapply(pred_files, 
                     function(x) gsub(".*/(.*)\\.tif", "\\1", x) %in% pred_vars) %>%
  pred_files[.]

# Confirm that all predictors are present
dummy <- sapply(pred_vars, function(x){
  times_present <- grepl(paste0("*./", x, "[\\./].*"),
                         c(ecodes_vars, pred_files)) %>%
    sum()
  if(times_present == 1){
    cat(x, "ok.\n")
  } else {
      warning(x, "is not present or present multiple times")
    }
  })
rm(dummy)

# get list of EcoDes-DK15 tiles (using dir here to speed up listing)
ecodes_tiles <- ecodes_vars[grepl("dtm", ecodes_vars)] %>%
  gsub("(.*/).*", "\\1", .) %>%
  gsub("/", "\\\\", .) %>%
  paste0("dir /b /s ", ., "*.tif") %>%
  shell(intern = T) %>%
  gsub("\\\\", "/", .)
ecodes_tiles <- gsub(".*([0-9]{4}_[0-9]{3}).tif", "\\1", ecodes_tiles)

# get mask folders /files
inland_mask <- ecodes_vars[grepl("inland_water", ecodes_vars)]
sea_mask <- ecodes_vars[grepl("sea", ecodes_vars)]
forest_mask_file <- "data/basemap_forests/forest_mask.tif"

# Write function to carry out predictions for one tile
project_tile <- function(tile_id, 
                         ecodes_preds,
                         pred_files,
                         inland_mask,
                         sea_mask,
                         forest_mask_file,
                         model_fit){
  # Status for debugging
  cat(tile_id, "\n")
  
  # Load ecodes rasters as stack
  tile_stack <- try(rast(paste0(ecodes_preds, 
                           gsub(".*/(.*)/$", "\\1", ecodes_preds),
                           "_", tile_id, ".tif")),
                    silent = T)
  # If that did not work we are dealing with a corner tile where the extent
  # differes between variables. We need to crop to common extent
  if(class(tile_stack) == "try-error"){
    # find smallest extent
    min_extent <- sapply(paste0(ecodes_preds, 
                       gsub(".*/(.*)/$", "\\1", ecodes_preds),
                       "_", tile_id, ".tif"),
           function(x) ncell(rast(x))) %>% which.min()
    min_extent <- ext(rast(paste0(ecodes_preds, 
                                     gsub(".*/(.*)/$", "\\1", ecodes_preds),
                                     "_", tile_id, ".tif")[min_extent]))
    tile_stack <- sapply(paste0(ecodes_preds, 
                                gsub(".*/(.*)/$", "\\1", ecodes_preds),
                                "_", tile_id, ".tif"),
                         function(x){
                           crop(rast(x), min_extent)
                         }) %>% Reduce(c, .)
    
  }
  names(tile_stack) <- gsub("_[0-9]{4}_[0-9]{3}", "", names(tile_stack))
  
  # Get tile boundaries
  tile_bounds <- ext(tile_stack)
  
  # Load an crop other predictors
  tile_stack <- sapply(pred_files, 
                        function(pred_file){
                          pred_raster <- rast(pred_file)
                          pred_raster <- crop(pred_raster, tile_bounds)
                          return(pred_raster)
                        }) %>%
    rast() %>%
    setNames(gsub(".*/(.*)\\..*", "\\1", pred_files)) %>%
    c(tile_stack, .)
  
  # Mask water from stack
  inland_mask_rast <- rast(paste0(inland_mask, "/inland_water_mask_", tile_id, ".tif"))
  if(ext(tile_stack) != ext(inland_mask_rast)) inland_mask_rast <- crop(inland_mask_rast, tile_bounds)
  tile_stack <- mask(tile_stack, inland_mask_rast)
  
  sea_mask_rast <- rast(paste0(sea_mask, "/sea_mask_", tile_id, ".tif"))
  if(ext(tile_stack) != ext(sea_mask_rast)) sea_mask_rast <- crop(sea_mask_rast, tile_bounds)
  tile_stack <- mask(tile_stack, sea_mask_rast)
  
  # Mask non forest from raster
  forest_mask <- rast(forest_mask_file)
  forest_mask_cropped <- crop(forest_mask, tile_bounds)
  tile_stack <- mask(tile_stack, forest_mask_cropped) 

  # project raster
  project_raster <- predict(as(tile_stack, "Raster"), model_fit)
  crs(project_raster) <- crs(tile_stack)
  
  # write out raster
  writeRaster(project_raster, 
               paste0("data/projections/gbm_biowide/forest_quality_" , tile_id, ".tif"),
               overwrite = T)
  
  # Return nothing
  return(NULL)
}
# project_tile(ecodes_tiles[12], ecodes_vars_sub, inland_mask, sea_mask, gbm_fit)

# Prepare parallel environmemt
cl <- makeCluster(46)
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(terra))
clusterEvalQ(cl, library(caret))
clusterEvalQ(cl, library(gbm))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(rgdal))
# Run projections
pblapply(ecodes_tiles,
         project_tile,
         ecodes_preds = ecodes_preds,
         pred_files = pred_files,
         inland_mask = inland_mask,
         sea_mask = sea_mask,
         forest_mask_file = forest_mask_file,
         model_fit = model_fit,
         cl = cl)

# Stop cluster
stopCluster(cl)

# Generate VRT file
setwd("data/projections/gbm_biowide/")
gdal_utils("buildvrt",
           source = list.files(".tif$"),
           destination = "forest_quality_gbm_biowide.vrt")
