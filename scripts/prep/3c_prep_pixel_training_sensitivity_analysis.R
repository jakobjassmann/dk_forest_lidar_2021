# DK Forest LiDAR - Data preparation script to generate training data for the
# sample size sensitivity analysis

# This is a modified version  of "3a_prep_pixel_training.R" developed in 
# response to reviewer 2's comments.

# Jakob Assmann jakob.assmann@uzh.ch 1 June 2024

# 1) Housekeeping ----

# Dependencies
library(terra)
library(raster)
library(parallel)
library(pbapply)
library(sf)
library(tidyverse)
library(purrr)
library(readxl)
library(factoextra)
library(ggplot2)
library(exactextractr)

# Set target number of samples
# (adding an extra buffer of 30% more to allow for removal of NAs due to water)
n_samples_list <- c(10000, 45000, 60000, 75000) * 1.3

for(n_samples in n_samples_list){
  cat(n_samples, "\n")
  
  # Load forest polygons
  load("data/training_data/polygon_geometries/high_quality_polys.Rda")
  load("data/training_data/polygon_geometries/low_quality_polys.Rda")
  
  ## 2) Generate pixel samples from forest polygons ----
  
  # Load EcoDes raster to determine pixels overlapping with polygons
  dtm_10m <- rast("F:/JakobAssmann/EcoDes-DK15_v1.1.0/dtm_10m/dtm_10m.vrt")
  
  # Setting the seed again for the sampling just in case the script crashed
  set.seed(414)
  
  # Retrieve coordinates and cell numbers of all pixels that touch the polygons
  high_quality_pixels <- exact_extract(dtm_10m, 
                                       high_quality, 
                                       include_xy = T, 
                                       include_cell = T)
  low_quality_pixels <- exact_extract(dtm_10m, 
                                      low_quality, 
                                      include_xy = T, 
                                      include_cell = T)
  
  # Add the polygon source to each of the pixels
  high_quality_pixels <- 1:nrow(high_quality) %>%
    pblapply(function(x){
      pixels <- high_quality_pixels[[x]]
      pixels$polygon_source <- high_quality[x,]$polygon_source
      return(pixels)
    })
  low_quality_pixels <- 1:nrow(low_quality) %>%
    pblapply(function(x){
      pixels <- low_quality_pixels[[x]]
      pixels$polygon_source <- low_quality[x,]$polygon_source
      return(pixels)
    })
  
  # First sample one pixel from each polygon
  high_quality_sample_ids <- high_quality_pixels %>%
    pbsapply(function(x) sample(x$cell, size = 1)) 
  low_quality_sample_ids <- low_quality_pixels %>%
    pbsapply(function(x) sample(x$cell, size = 1)) 
  
  # Next bind samples into one big data frame:
  high_quality_pixels <- bind_rows(high_quality_pixels)
  low_quality_pixels <- bind_rows(low_quality_pixels)
  
  # Extract and remove pixels already sampled
  high_quality_sample <- filter(high_quality_pixels, 
                                cell %in% high_quality_sample_ids)
  high_quality_pixels <- filter(high_quality_pixels, 
                                !(cell %in% high_quality_sample_ids))
  low_quality_sample <- filter(low_quality_pixels, 
                               cell %in% low_quality_sample_ids)
  low_quality_pixels <- filter(low_quality_pixels, 
                               !(cell %in% low_quality_sample_ids))
  
  # Finally sample sample reamining pixels to make up a sample of 30k pixels
  # each:
  high_quality_sample <- sample_n(high_quality_pixels,
                                  n_samples - nrow(high_quality_sample)) %>%
    bind_rows(high_quality_sample, .)
  low_quality_sample <- sample_n(low_quality_pixels,
                                 n_samples - nrow(low_quality_sample)) %>%
    bind_rows(low_quality_sample, .)
  
  # Convert dataframes for sf objects and clean up
  high_quality_sample <- st_as_sf(high_quality_sample,
                                  coords = c("x", "y"),
                                  crs = st_crs(dtm_10m)) %>%
    select(cell, polygon_source, geometry)
  low_quality_sample <- st_as_sf(low_quality_sample,
                                 coords = c("x", "y"),
                                 crs = st_crs(dtm_10m)) %>%
    select(cell, polygon_source, geometry)
  
  # Add forest quality and sample id columns
  high_quality_sample <- high_quality_sample %>%
    mutate(forest_value = "high",
           sample_id = paste0("high_", 1:nrow(high_quality_sample)))
  low_quality_sample <- low_quality_sample %>%
    mutate(forest_value = "low",
           sample_id = paste0("low_", 1:nrow(low_quality_sample)))
  
  
  # Check stats
  high_quality_sample %>%
    st_drop_geometry() %>%
    group_by(polygon_source) %>% 
    tally()
  low_quality_sample %>% 
    st_drop_geometry() %>%
    group_by(polygon_source) %>% 
    tally()
  
  # Save / load
  save(high_quality_sample,
       file = paste0("data/training_data/pixel_geometries/high_quality_pixel_sample_", n_samples, ".Rda"))
  save(low_quality_sample, 
       file = paste0("data/training_data/pixel_geometries/low_quality_pixel_sample_", n_samples, ".Rda"))

  # bind into one sf object 
  combined_sample_coords <- rbind(high_quality_sample,
                                  low_quality_sample)
  
  # Save / load
  save(combined_sample_coords, 
       file = paste0("data/training_data/pixel_geometries/combined_pixel_sample_", n_samples, ".Rda"))
  # load("data/training_data/pixel_geometries/combined_pixel_sample.Rda")
  
  # Convert to vect (terra) for fast extraction
  combined_sample_coords <- combined_sample_coords %>%
    vect()
  
  ## 4) Extract predictor variables for training locations ----
  
  ## EcoDes-DK15 v1.1.0 descriptors
  
  # Load list of Ecodes-DK variables
  ecodes_vrt <- shell("dir /b /s F:\\JakobAssmann\\EcoDes-DK15_v1.1.0\\*.vrt",
                      intern = T) %>%
    gsub("\\\\", "/", .)
  
  # Remove all variables not relevant
  ecodes_vrt <- ecodes_vrt[!grepl("date_stamp",ecodes_vrt)]
  ecodes_vrt <- ecodes_vrt[!grepl("point_count",ecodes_vrt)]
  ecodes_vrt <- ecodes_vrt[!grepl("point_source",ecodes_vrt)]
  ecodes_vrt <- ecodes_vrt[!grepl("building",ecodes_vrt)]
  
  # Write helper function to extract samples
  extract_fun <- function(vrt_file){
    # Read in raster
    vrt_rast <- terra::rast(vrt_file)
    
    # Extract 
    cell_values <- terra::extract(vrt_rast, combined_sample_coords)[,2]
    
    # Prepare extracted values for return as data.frame
    extractions <- data.frame(
      sample_id = combined_sample_coords$sample_id,
      forest_value = combined_sample_coords$forest_value,
      cell_values = cell_values)
    
    # Update colum column name with name of vrt files
    names(extractions)[3] <- gsub(".*/(.*).vrt", "\\1", vrt_file)
    return(extractions)
  }
  
  # Prepare cluster
  cl <- makeCluster(11)
  clusterEvalQ(cl, library(terra))
  
  # Export coordinates
  combined_sample_coords_wrapped <- wrap(combined_sample_coords)
  clusterExport(cl, varlist = "combined_sample_coords_wrapped")
  clusterEvalQ(cl, {
    combined_sample_coords <- vect(combined_sample_coords_wrapped)
    print(head(combined_sample_coords))
  })
  
  # Extract variables in parallel (took around 1h with 30 cores on d23510)
  combined_sample <- pblapply(ecodes_vrt,
                              extract_fun,
                              cl = cl) 
  
  save(combined_sample, 
       file = paste0("data/training_data/ecodes_pixel_sample_temp_", n_samples, ".Rda"))
  
  # Stop cluster
  stopCluster(cl)
  rm(cl)
  
  # Combine into one single dataframe and merge with sf object to add coordinates 
  # to sample
  pixel_training_data_raw <- combined_sample %>% 
    map(function(x) dplyr::select(x, -forest_value)) %>%
    reduce(full_join, by = "sample_id") %>% 
    full_join(rbind(high_quality_sample, 
                    low_quality_sample), 
              ., by = "sample_id")
  
  # confirm order is the same as in the combined_sample_coords vect object
  nrow(combined_sample_coords)
  sum(combined_sample_coords$sample_id == pixel_training_data_raw$sample_id)
  head(combined_sample_coords$sample_id == pixel_training_data_raw$sample_id)
  tail(combined_sample_coords$sample_id == pixel_training_data_raw$sample_id)
  which(!(combined_sample_coords$sample_id == pixel_training_data_raw$sample_id))
  
  ## Bjerreskov et al. 2021 Forest type (coniferous / decidious)
  
  # Load rasters
  treetype_bjer_dec <- rast("data/predictor_data/treetype/treetype_bjer_dec.tif")
  treetype_bjer_con <- rast("data/predictor_data/treetype/treetype_bjer_con.tif")
  
  # Extract treetype
  pixel_training_data_raw$treetype_bjer_dec <- terra::extract(treetype_bjer_dec, combined_sample_coords)[,2]
  pixel_training_data_raw$treetype_bjer_con <- terra::extract(treetype_bjer_con, combined_sample_coords)[,2]
  
  # Unload rasters
  rm(treetype_bjer_dec)
  rm(treetype_bjer_con)
  
  ## Focal variables
  
  # Get list of rasters
  focal_vars <- list.files("data/predictor_data/focal_variables/", "tif", full.names = T)
  
  # Prepare cluster
  cl <- makeCluster(16)
  clusterEvalQ(cl, library(terra))
  
  # Export coordinates
  combined_sample_coords_wrapped <- wrap(combined_sample_coords)
  clusterExport(cl, varlist = "combined_sample_coords_wrapped")
  clusterEvalQ(cl, {
    combined_sample_coords <- vect(combined_sample_coords_wrapped)
    print(head(combined_sample_coords))
  })
  
  # Extract variables in parallel (took around 34 s on d23510)
  focal_vars <- pblapply(focal_vars,
                         function(rast_file){
                           focal_var <- rast(rast_file)
                           cell_values <- terra::extract(focal_var, combined_sample_coords)[,2]
                           extractions <- data.frame(
                             sample_id = combined_sample_coords$sample_id,
                             cell_values = cell_values)
                           names(extractions)[2] <-  gsub(".*/(.*).tif", "\\1", rast_file) 
                           return(extractions)
                         },
                         cl = cl) %>%
    reduce(full_join, by = "sample_id") 
  
  # Stop cluster
  stopCluster(cl)
  rm(cl)
  
  # Merge outputs
  pixel_training_data_raw <- full_join(pixel_training_data_raw, focal_vars, by = "sample_id")
  
  ## Near surface groundwater(summer)
  
  # Load raster
  ns_groundwater_summer_utm32_10m <- rast("data/predictor_data/terraennaert_grundvand_10m/ns_groundwater_summer_utm32_10m.tif")
  
  # Extract values
  pixel_training_data_raw$ns_groundwater_summer_utm32_10m <- terra::extract(ns_groundwater_summer_utm32_10m, combined_sample_coords)[,2]
  
  # Remove raster
  rm(ns_groundwater_summer_utm32_10m)
  
  ## Soil variables from Derek / SustainScapes 
  
  # Load raster
  Clay_utm32_10m <- rast("data/predictor_data/soil_layers/Clay_utm32_10m.tif")
  Sand_utm32_10m <- rast("data/predictor_data/soil_layers/Sand_utm32_10m.tif")
  Soc_utm32_10m <- rast("data/predictor_data/soil_layers/Soc_utm32_10m.tif")
  
  # Extract values
  pixel_training_data_raw$Clay_utm32_10m <- terra::extract(Clay_utm32_10m, combined_sample_coords)[,2]
  pixel_training_data_raw$Sand_utm32_10m <- terra::extract(Sand_utm32_10m, combined_sample_coords)[,2]
  pixel_training_data_raw$Soc_utm32_10m <- terra::extract(Soc_utm32_10m, combined_sample_coords)[,2]
  
  # Remove rasters
  rm(list = c("Clay_utm32_10m", "Sand_utm32_10m", "Soc_utm32_10m"))
  
  ## Foliage height diversity
  
  # Load raster
  foliage_height_diversity <- rast("data/predictor_data/foliage_height_diversity/foliage_height_diversity.tif")
  
  # Extract values
  pixel_training_data_raw$foliage_height_diversity <- terra::extract(foliage_height_diversity, combined_sample_coords)[,2]
  
  # Remove raster
  rm(foliage_height_diversity)
  
  # Save intermediate backup
  save(pixel_training_data_raw, 
       file = paste0("data/training_data/pixel_training_", n_samples, ".Rda"))
  # load("data/training_data/pixel_training.Rda")
  
  ## 5) Add stratification
  
  ## BIOWIDE stratification
  
  # Load geometries
  biowide_strat <- read_sf("data/stratification/biowide_georegions/biowide_zones.shp")
  
  # Extract regions
  pixel_training_data_raw$biowide_region <-
    pixel_training_data_raw %>%
    st_as_sf() %>%
    st_intersects(biowide_strat) %>%
    split(., seq_along(.)) %>%
    lapply(function(x){
      x <- unlist(x)
      # Check if there is no intersection, if so set to dummy region outside region range (will produce NAs in classification)
      if(length(x) == 0) x <- list(c(nrow(biowide_strat) + 1))
      return(x)}) %>% 
    unlist() %>%
    biowide_strat$region[.]
  
  ## Derek's stratification
  dereks_strat <- rast("data/stratification/derek_stratification/Results_2clim_5soil.tif")
  
  # Extract data
  pixel_training_data_raw$dereks_stratification <- terra::extract(dereks_strat, combined_sample_coords)[,2]
  
  
  # Save final version
  pixel_training_data <- st_as_sf(pixel_training_data_raw)
  pixel_training_data <- relocate(pixel_training_data, 
                                  forest_value,
                                  sample_id,
                                  polygon_source, 
                                  biowide_region,
                                  dereks_stratification)
  save(pixel_training_data, file = paste0("data/training_data/pixel_training_", n_samples, ".Rda"))
  
  rm(list = ls())
}  

