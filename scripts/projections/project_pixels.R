# DK Forest LiDAR by pixel projections across the EcoDes dataset
# Jakob Assmann j.assmann@bio.au.dk 9 July 2021

## 1) Prep environment ----

# Dependencies
library(caret)
library(gbm)
library(raster)
library(terra)
library(parallel)
library(pbapply)
library(sf)
library(rgdal)

# Load models
load("data/models/final_gbm_model_pixel_biowide.Rda")
load("data/models/final_ranger_model_pixel_biowide.Rda")

## 2) Function definitions ----
# Write function to carry out predictions for one tile
project_tile <- function(tile_id, 
                         ecodes_preds,
                         pred_files,
                         inland_mask,
                         sea_mask,
                         forest_mask_file,
                         model_fit,
                         model_name){
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
               paste0("data/projections/", model_name, 
                      "/forest_quality_" , tile_id, ".tif"),
               overwrite = T)
  
  # Return nothing
  return(NULL)
}
# Function test line:
# project_tile(ecodes_tiles[12], ecodes_vars_sub, inland_mask, sea_mask, gbm_fit)

# Function to project all tiles for a model
project_model <- function(model_fit,
                          model_name,
                          ncores){
  # Status
  cat("Projections for", model_name, "\n\n")
  
  # Check whether folder for projections exists, prompt user in case it does
  # and delete if requested. 
  if(dir.exists(paste0("data/projections/", model_name, "/"))){
    delete <- readline(paste0("Folder exists: 'data/projections/", model_name, 
                              "/' - Delete folder? [y/n]: "))
    if(delete == "y") { 
      unlink(paste0("Folder exists: data/projections/"), 
             recursive=TRUE)
    } else if(delete == "n"){
      cat("\nOkay, stopping projections.")
      return(1)
    } else{
      cat("\nInvalid response, stopping projections.")
      return(1)
    }
  }
  # Create dir
  dir.create(paste0("data/projections/", model_name, "/"))
  
  # Prepare environment for model projections
  cat("Gathering EcoDes-DK15 predictor layers ...\n")
  
  # Get list of predictor variables 
  if(model_fit$method == "gbm"){
    pred_vars <- summary(model_fit, plotit = F)$var
  } else if(model_fit$method == "ranger") {
    pred_vars <- row.names(varImp(model_fit)$importance)
  } else {
    warning("Model Type:", model_fit$method, "not recognised.",
            "\n Adapt project_model() to method to allow for predictor",
            "variable extraction.")
  }   
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
  cat("Gathering additional predictor layers ...\n")
  pred_files <- list.files("data/predictor_data/", "tif", recursive = T, 
                           full.names = T)
  pred_files <- sapply(pred_files, 
                       function(x) gsub(".*/(.*)\\.tif", "\\1", x) %in% pred_vars) %>%
    pred_files[.]
  
  # get list of EcoDes-DK15 tiles (using dir here to speed up listing)
  ecodes_tiles <- ecodes_vars[grepl("dtm", ecodes_vars)] %>%
    gsub("(.*/).*", "\\1", .) %>%
    gsub("/", "\\\\", .) %>%
    paste0("dir /b /s ", ., "*.tif") %>%
    shell(intern = T) %>%
    gsub("\\\\", "/", .)
  ecodes_tiles <- gsub(".*([0-9]{4}_[0-9]{3}).tif", "\\1", ecodes_tiles)

  # Confirm that all predictors are present
  cat("\nConfirming  all predictor layers are present:\n\n")
  dummy <- sapply(pred_vars, function(x){
    times_present <- grepl(paste0("*./", x, "[\\./].*"),
                           c(ecodes_vars, pred_files)) %>%
      sum()
    if(times_present == 1){
      cat("\t", x, "\tok.\n")
    } else {
      warning(x, "is not present or present multiple times")
      return(1)
    }
  })
  rm(dummy)
  
  # get mask folders /files
  cat("\nGathering masks ...\n")
  inland_mask <- ecodes_vars[grepl("inland_water", ecodes_vars)]
  sea_mask <- ecodes_vars[grepl("sea", ecodes_vars)]
  forest_mask_file <- "data/basemap_forests/forest_mask.tif"
  
  # Prepare parallel environmemt
  cat("Preparing parallel environment ...\n")
  cl <- makeCluster(ncores)
  clusterEvalQ(cl, library(raster))
  clusterEvalQ(cl, library(terra))
  clusterEvalQ(cl, library(caret))
  clusterEvalQ(cl, library(gbm))
  clusterEvalQ(cl, library(dplyr))
  clusterEvalQ(cl, library(rgdal))
  
  # Run projections
  cat("\nStarting projections...\n")
  pblapply(ecodes_tiles,
           project_tile,
           ecodes_preds = ecodes_preds,
           pred_files = pred_files,
           inland_mask = inland_mask,
           sea_mask = sea_mask,
           forest_mask_file = forest_mask_file,
           model_fit = model_fit,
           model_name = model_name,
           cl = cl)
  
  # Stop cluster
  cat("Stopping cluster ...\n")
  stopCluster(cl)

  current_wd <- getwd()
  setwd(paste0("data/projections/", model_name, "/"))
  
  # Generate VRT file
  cat("Generating VRT file ... \n")
  gdal_utils("buildvrt",
             source = list.files(pattern = "tif$"),
             destination = paste0("forest_quality_", model_name, ".vrt"))
  
  # Generate Cloud Optimised GeoTiff from VRT (and project to EPSG:4326)
  cat("Generating Cloud Optimised GeoTiff ... \n")
  gdal_utils("warp",
             source = paste0("forest_quality_", model_name, ".vrt"),
             destination = paste0("forest_quality_", model_name, "_cog_epsg4326.tif.vrt"),
             options = c(
               "-t_srs", "EPSG:4326",
               "-of", "COG",
               "-co", "RESAMPLING=NEAREST",
               "-co", "TILING_SCHEME=GoogleMapsCompatible",
               "-co", "COMPRESS=DEFLATE",
               "-co", "NUM_THREADS=46"
             ))
  setwd(current_wd)
  
  # Status
  cat("Done.\n\n")
  
  # Return
  return(NULL)
}

## 3) Execute projections
project_model(gbm_fit, "gbm_biowide", 46)
project_model(rf_fit, "ranger_biowide", 46)
## End of file
