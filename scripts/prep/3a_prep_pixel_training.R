# DK Forest LiDAR - Data preparation script to generate training data
# Jakob Assmann j.assmann@bio.au.dk 4 February 2022

# 1) Housekeeping ----

# Dependencies
library(terra)
library(raster)
library(parallel)
library(pbapply)
library(sf)
library(tidyverse)
library(purrr)
library(dplyr)
library(readxl)
library(factoextra)
library(ggplot2)

# Set seed for pseudo random numbers
set.seed(1168)

## 2) Prepare forest training polygons ----

# Load geometries for high quality forests
high_quality <- list.files("data/response_data/high_quality_forests/", 
                           ".shp", 
                           recursive = T,
                           full.names = T) %>%
  # Load files and assign source column
  map(function(shp_file){
    polygons <- read_sf(shp_file)
    # Filter old growth forests if needed (aftaler om natur)
    if(sum("tilskudsor" %in% names(polygons)) > 0){
      polygons <- filter(polygons, tilskudsor == "Privat urørt skov")
    }
    polygons <- select(polygons, !everything())
    # make geometries valid remove have a 10 m pixel diagonal and any overlap
    polygons <- polygons %>% 
      st_make_valid() %>%
      st_buffer(-(sqrt(10^2+10^2)/2)) %>%
      st_difference()
    # Add polygons source 
    polygons$polygon_source <- gsub(".*/(.*)\\.shp", "\\1", shp_file)
    return(polygons)
    }) %>%
  bind_rows() %>%
  mutate(polygon_source = case_when(
    polygon_source == "skov_kortlaegning_2016_2018" ~ "p15",
    polygon_source == "p25_offentligareal" ~ "p25",
    polygon_source == "aftale_natur_tinglyst" ~ "private_old_growth",
  )) 

## The next section is to remove overlapping geometries from the high quality data

# Re arrange order of dataframe so that p25 polygons are on top
high_quality <- high_quality[order(match(high_quality$polygon_source, c("p25", "private_old_growth", "p15"))),]

# add arbitary id
high_quality$id <- paste0(high_quality$polygon_source, "_", 1:nrow(high_quality))
# Identify polygons that overlap and which other polygon they overlap with
high_quality$overlaps_with <- high_quality %>%
  st_overlaps()
high_quality <- high_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
# Check out stats
high_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()

# Clean overlap starting with the p25 polygons as they are first up
for(poly_id in high_quality$id[high_quality$overlaps]){
  cat("Cleaning up for", poly_id, "\n")
  # Get polygon geometry
  poly_geo <- high_quality %>% 
    filter(id == poly_id)
  overlaps_with <- high_quality %>% 
    filter(id == poly_id) %>%
    pull(overlaps_with) %>%
    unlist() 
  # Remove from geometry from all other polygons that it overlaps with
  high_quality[overlaps_with,]$geometry <- overlaps_with %>% 
    slice(high_quality, .) %>%
    split(., 1:nrow(.)) %>%
    map(function(x) st_difference(st_geometry(x), poly_geo)) %>%
    do.call(c, .)
  # Clean up environment
  rm(poly_geo)
  rm(overlaps_with)
}

# Check results
high_quality$overlaps_with <- high_quality %>%
  st_overlaps()
high_quality <- high_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
high_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()

# Now for some reason some polygons with a line string overlap will still remain
# Here we do a short cut by shrinking those by 0.1 m. 
high_quality[high_quality$overlaps,]$geometry <- high_quality %>% 
  filter(overlaps) %>%
  split(., 1:nrow(.)) %>%
  map(function(x) st_buffer(x, -0.1)) %>%
  bind_rows() %>% 
  st_geometry()

# Confirm that worked
high_quality$overlaps_with <- high_quality %>%
  st_overlaps()
high_quality <- high_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
high_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()

# Remove excess colums
high_quality <- high_quality %>%
  select(polygon_source, geometry)

## Load and prep geometries for low quality forests
# Plantations geometries and meta data
plantations <- read_sf("data/response_data/low_quality_forests/NST_plantations/LitraPolygoner_region/LitraPolygoner_region.shp")
plantations_meta <- read_excel("data/response_data/low_quality_forests/NST_plantations/NST  2019 08012019 ber 16012020 til bios_au.xlsx") 
# Helper function to classify age bins (given as decadal mid-points, e.g., 5, 15, 25, 35 etc.)
sort_into_age_bins <- function(Aldersklasse){
  case_when(Aldersklasse > 0 & Aldersklasse <= 5 ~ "0_to_10",  # Filters decadal mid-point 5
            Aldersklasse > 5 & Aldersklasse <= 25 ~ "10_to_30", # Filters decadal mid-points 15 and 25
            Aldersklasse > 25 & Aldersklasse <= 45 ~ "30_to_50", # Filters decadal mid-points 35 and 45
            Aldersklasse > 45 & Aldersklasse <= 65 ~ "50_to_70", # Filters decadal mid-poinst 55 and 65
            Aldersklasse > 65 & Aldersklasse <= 95 ~ "70_to_100", # Filters decadal mid-points 75, 85 and 95
            TRUE ~ "NA") # Everything else is set to NA
}
# Data cleaning
plantations_meta <- plantations_meta %>%
  select(Ident, Aldersklasse, `ANV 3`, `Allerede urørt`, Status) %>%
  filter(`ANV 3` == 0) %>% # Keep only ANV 3 values of “0” (remove everything else)
  filter(`Allerede urørt` != "Urørt") %>% # Throw out all rows with label “Urørt"
  filter(Status == "G") %>% # Keep only current plantations (status = “G”) 
  na.omit() %>%
  mutate(age_bin = sort_into_age_bins(Aldersklasse)) %>%
  filter(age_bin != "NA") %>% # Remove all age classes not in the above categories
  group_by(age_bin) %>%
  sample_n(1000)
# Filter geometries
plantations <- plantations %>%
  st_make_valid() %>%
  filter(UNIKID %in% plantations_meta$Ident)
# Set source column
plantations$polygon_source <- "NST_plantations"
# Remove all unwanted columns
plantations <- select(plantations, geometry, polygon_source)
# Buffer geometries by half a pixel diagonal
plantations <- st_buffer(plantations, -(sqrt(10^2+10^2)/2))

# Load ikke p25 geometries
ikke_p25 <- read_sf("data/response_data/low_quality_forests/ikke_p25/ikkeP25_skov.shp") 
ikke_p25 <- ikke_p25 %>%
  st_make_valid() 
ikke_p25$polygon_source <- "ikke_p25"
# Remove all unwanted columns
ikke_p25 <- select(ikke_p25, geometry, polygon_source)

# Buffer geometries by half a pixel diagonal
ikke_p25 <- st_buffer(ikke_p25, -(sqrt(10^2+10^2)/2))

# Remove overlap between high quality and low quality polygons,
# as well as overlap between low quality polygons

# first we append the low quality polygons to the high quality data
# we start with the ikke_p25 forest to prioritise those over the plantations
low_quality <- bind_rows(high_quality,
                         ikke_p25,
                         plantations)

# add a unique polygon id
low_quality$id <- paste0(low_quality$polygon_source, "_", 1:nrow(low_quality))

# Identify polygons that overlap and which other polygon they overlap with
low_quality$overlaps_with <- low_quality %>%
  st_overlaps()
low_quality <- low_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
# Check out stats
low_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()
# # A tibble: 10 x 2
# n_overlaps     n
# <int> <int>
#   1          0 17341
# 2          1  2043
# 3          2   314
# 4          3    75
# 5          4    24
# 6          5     8
# 7          6     5
# 8          7     1
# 9          8     3
# 10         11     1

# Remove overlap starting from top to bottom working only on low quality polys
for(poly_id in low_quality$id[low_quality$overlaps]){
  cat("Cleaning up for", poly_id, "\n")
  # Get polygon geometry
  poly_geo <- low_quality %>% 
    filter(id == poly_id)
  overlaps_with <- low_quality %>% 
    filter(id == poly_id) %>%
    pull(overlaps_with) %>%
    unlist() 
  # Remove all polygons high_quality polygons
  overlaps_with <- overlaps_with[!(low_quality$polygon_source[overlaps_with] %in% c("p25", "private_old_growth", "p15"))]
  # Remove from geometry from all other polygons that it overlaps with
  if(length(overlaps_with) > 0) {
    low_quality[overlaps_with,]$geometry <- overlaps_with %>% 
    slice(low_quality, .) %>%
    split(., 1:nrow(.)) %>%
    map(function(x){
      difference <- st_difference(st_geometry(x), poly_geo)
      if(length(difference) == 0) difference <- st_geometry(st_buffer(poly_geo, -10^6))
      return(difference)
    }) %>%
    do.call(c, .)
  }
  # Clean up environment
  rm(poly_geo)
  rm(overlaps_with)
}

# Remove empty geometries
low_quality <- filter(low_quality, !st_is_empty(low_quality))

# Check whether this worked
low_quality$overlaps_with <- low_quality %>%
  st_overlaps()
low_quality <- low_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
# Check out stats
low_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()
# # A tibble: 10 x 2
# n_overlaps     n
# <int> <int>
#   1          0 17619
# 2          1  1805
# 3          2   211
# 4          3    49
# 5          4     9
# 6          5     7
# 7          6     2
# 8          7     1
# 9          8     1
# 10         10     1
# -> Not bad it definitely improved it. 

# The remaining overlaps are likely line overlaps. 
# Let's see whether we can address the issue by shrinking all affected low
# quality polygons by -0.1 m
low_quality <- low_quality %>% 
  filter(polygon_source %in% c("ikke_p25", "NST_plantations" )) %>%
  filter(overlaps) %>%
  split(., 1:nrow(.)) %>%
  map(function(x) st_buffer(x, -0.1)) %>%
  bind_rows() %>% 
  bind_rows(filter(low_quality, !(polygon_source %in% c("ikke_p25", "NST_plantations" ))),
            filter(low_quality, (polygon_source %in% c("ikke_p25", "NST_plantations" ) & !overlaps)),
            .)

# Check whether this worked
low_quality$overlaps_with <- low_quality %>%
  st_overlaps()
low_quality <- low_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
# Check out stats
low_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()
# # A tibble: 2 x 2
# n_overlaps     n
# <int> <int>
#   1          0 19703
# 2          1     2
# Only two more forests!
# Check out the remaining forests
test <- low_quality[low_quality$overlaps,]
plot(st_geometry(test[1,]))
plot(st_geometry(test[2,]), add = T, col = "red")
# The forests are identical, remove the ikke_p25 forest and finalise low_quality 
# set by removing the high quality forests
low_quality <- filter(low_quality, id != (low_quality[low_quality$overlaps,] %>% 
                               filter(polygon_source == "ikke_p25") %>%
                               pull(id))) %>%
  filter(polygon_source %in% c("ikke_p25", "NST_plantations"))
# Remove surplus columns
low_quality <- select(low_quality, polygon_source, geometry)


# Final check to make sure there is no overlap
st_overlaps(high_quality) %>% sapply(length) %>% sum()
st_overlaps(low_quality) %>% sapply(length) %>% sum()
st_overlaps(high_quality, low_quality) %>% sapply(length) %>% sum()

# Save / load geometries
save(high_quality, file = "data/training_data/polygon_geometries/high_quality_polys.Rda")
save(low_quality, file = "data/training_data/polygon_geometries/low_quality_polys.Rda")
# load("data/training_data/polygon_geometries/high_quality_polys.Rda")
# load("data/training_data/polygon_geometries/low_quality_polys.Rda")

## 3) Generate pixel samples from forest polygons ----

# First, if possible, sample every polygon once at random within 10 m distance 
# from edge 
high_quality_sample_coords <- high_quality %>%
  split(., 1:nrow(.)) %>%
  pblapply(st_sample, size = 1) %>%
  do.call(c, .) %>%
  st_sf() %>%
  mutate(forest_value = "high") %>% 
  mutate(sample_id = paste0(forest_value, "_", 1:n())) 
low_quality_sample_coords <- low_quality %>%
  split(., 1:nrow(.)) %>%
  pblapply(st_sample, size = 1) %>%
  do.call(c, .) %>%
  st_sf() %>%
  mutate(forest_value = "low") %>% 
  mutate(sample_id = paste0(forest_value, "_", 1:n())) 

# Next sample across all polygons at random to make up 30k samples
high_quality_sample_coords_rand <- st_sample(
  high_quality, 
  30000 - nrow(high_quality_sample_coords)) %>% 
  st_sf() %>%
  mutate(forest_value = "high") %>% 
  mutate(sample_id = paste0(forest_value, "_", nrow(high_quality_sample_coords) + 1:n())) 
low_quality_sample_coords_rand <- st_sample(
  low_quality, 
  30000 - nrow(low_quality_sample_coords)) %>% 
  st_sf() %>%
  mutate(forest_value = "low") %>% 
  mutate(sample_id = paste0(forest_value, "_", nrow(high_quality_sample_coords) + 1:n())) 

# Load EcoDes raster to allow for checking of duplicate pixels
dtm_10m <- rast("F:/JakobAssmann/EcoDes-DK_v1.1.0/dtm_10m/dtm_10m.vrt")

# Define helper function to find duplicate pixels
find_duplicates <- function(first_sample, second_sample = NULL){
  # Combine samples if needed
  if(!is.null(second_sample)){
    both_samples <- bind_rows(first_sample, second_sample) %>% vect()
  } else {
    # Otherwise assign and duplicate first sample for convenience subsetting later
    both_samples <- vect(first_sample)
    second_sample <- first_sample
  }
  # Extract pixels from dtm raster to get cell ID
  test_extract <- terra::extract(dtm_10m, both_samples, cells = T)
  # Check whether any of the cell ids sampled are present more than twice
  duplicates <-  
    pbsapply(test_extract$cell, function(x) sum(test_extract$cell %in% x) > 1)
  # Associated cell ids with sample_ids
  duplicates <- both_samples[test_extract$ID[duplicates],]$sample_id
  # Return subset of second sample containing only the duplicated ids
  return(
    filter(second_sample, 
           sample_id %in% duplicates))
}

# Keep searching for duplicates till none are left

# Dummy dataframe
high_quality_duplicates <- data.frame(1:0)
counter <- 1
while(nrow(high_quality_duplicates) > 0){
  # Status
  cat("Removing high quality sample duplicates - round: ", counter, "\n",
      "Duplicate samples to replace: ", nrow(high_quality_duplicates), "\n")
  # Use helper function to find duplicates
  high_quality_duplicates <- find_duplicates(
    high_quality_sample_coords,
    high_quality_sample_coords_rand)
  # Re-extract duplicates
  high_quality_sample_coords_rand <- st_sample(
    high_quality, size = nrow(high_quality_duplicates)) %>% 
    st_sf() %>%
    mutate(forest_value = "high") %>% 
    mutate(sample_id = high_quality_duplicates$sample_id) %>%
    # Remove old duplicates and add replacement coordinates
    bind_rows(filter(high_quality_sample_coords_rand,
                     !(sample_id %in% high_quality_duplicates$sample_id),
                     ), .)
  counter <- counter + 1
}
# Merge final samples
high_quality_sample_coords <- bind_rows(high_quality_sample_coords,
                                        high_quality_sample_coords_rand)

# Dummy dataframe
low_quality_duplicates <- data.frame(1:0)
counter <- 1
while(nrow(low_quality_duplicates) > 0){
  # Status
  cat("Removing low quality sample duplicates - round: ", counter, "\n",
      "Duplicate samples to replace: ", nrow(low_quality_duplicates), "\n")
  # Use helper function to find duplicates
  low_quality_duplicates <- find_duplicates(
    low_quality_sample_coords,
    low_quality_sample_coords_rand)
  # Re-extract duplicates
  low_quality_sample_coords_rand <- st_sample(
    low_quality, size = nrow(low_quality_duplicates)) %>% 
    st_sf() %>%
    mutate(forest_value = "low") %>% 
    mutate(sample_id = low_quality_duplicates$sample_id) %>%
    # Remove old duplicates and add replacement coordinates
    bind_rows(filter(low_quality_sample_coords_rand,
                     !(sample_id %in% low_quality_duplicates$sample_id),
    ), .)
  counter <- counter + 1
}
# Merge final samples
low_quality_sample_coords <- bind_rows(low_quality_sample_coords,
                                        low_quality_sample_coords_rand)


# Add polygon source
high_quality_sample_coords$polygon_source <- st_intersects(
  high_quality_sample_coords,
  high_quality) %>%
  sapply(function(x) high_quality$polygon_source[x[1]])
low_quality_sample_coords$polygon_source <- st_intersects(
  low_quality_sample_coords,
  low_quality) %>%
  sapply(function(x) low_quality$polygon_source[x[1]])

# Save / load
save(high_quality_sample_coords, file = "data/training_data/pixel_geometries/high_quality_pixel_sample.Rda")
save(low_quality_sample_coords, file = "data/training_data/pixel_geometries/low_quality_pixel_sample.Rda")
# load("data/training_data/pixel_geometries/high_quality_pixel_sample.Rda")
# load("data/training_data/pixel_geometries/low_quality_pixel_sample.Rda")

# transform geometries to raster ETRS89 / UTM32 
target_crs <- st_crs(raster("F:/JakobAssmann/EcoDes-DK_v1.1.0/dtm_10m/dtm_10m_6049_684.tif"))
high_quality_sample_coords <- high_quality_sample_coords %>%
  st_transform(target_crs)
low_quality_sample_coords <- low_quality_sample_coords %>%
  st_transform(target_crs)

# bind into one sf object 
combined_sample_coords <- rbind(high_quality_sample_coords,
                                low_quality_sample_coords)

# Save / load
save(combined_sample_coords, file = "data/training_data/pixel_geometries/combined_pixel_sample.Rda")
# load("data/training_data/pixel_geometries/combined_pixel_sample.Rda")

# Write as shp (not required for parallel processing with terra anymore)
# write_sf(combined_sample_coords, 
#          dsn = "data/training_data/pixel_geometries/combined_pixel_sample.shp")

# Convert to vect (terra) for fast extraction
combined_sample_coords <- combined_sample_coords %>%
  as_Spatial() %>%
  vect()

## 4) Extract predictor variables for training locations ----

## EcoDes-DK15 v1.1.0

# Load list of Ecodes-DK variables
ecodes_vrt <- shell("dir /b /s F:\\JakobAssmann\\EcoDes-DK_v1.1.0\\*.vrt",
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
    forest_value = combined_sample_coords$forest_value ,
    cell_values = cell_values)
  
  # Update colum column name with name of vrt files
  names(extractions)[3] <- gsub(".*/(.*).vrt", "\\1", vrt_file)
  return(extractions)
}

# Prepare cluster
cl <- makeCluster(30)
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

save(combined_sample, file = "data/training_data/ecodes_pixel_sample_temp.Rda")
# load("data/training_data/ecodes_pixel_sample_temp.Rda")

# Stop cluster
stopCluster(cl)
rm(cl)

# Combine into one single dataframe and merge with sf object to add coordinates 
# to sample
pixel_training_data_raw <- combined_sample %>% 
  map(function(x) dplyr::select(x, -forest_value)) %>%
  reduce(full_join, by = "sample_id") %>% 
  full_join(rbind(high_quality_sample_coords, 
                  low_quality_sample_coords), 
            ., by = "sample_id")

# confirm order is the same as in the combined_sample_coords vect object
nrow(combined_sample_coords)
sum(combined_sample_coords$sample_id == pixel_training_data_raw$sample_id)
head(combined_sample_coords$sample_id == pixel_training_data_raw$sample_id)
tail(combined_sample_coords$sample_id == pixel_training_data_raw$sample_id)
which(!(combined_sample_coords$sample_id == pixel_training_data_raw$sample_id))

# ## Alex and Jakob's forest type (coniferous vs. broadleaf)
# 
# # Load rasters
# forest_type_cloud <- rast("data/predictor_data/conif_vs_broadleaf/forest_type_cloud.tif")
# forest_type_con <- rast("data/predictor_data/conif_vs_broadleaf/forest_type_con.tif")
# forest_type_dec <- rast("data/predictor_data/conif_vs_broadleaf/forest_type_dec.tif")
# bornholm_forest_type_cloud <- rast("data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_cloud.tif")
# bornholm_forest_type_con <- rast("data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_con.tif")
# bornholm_forest_type_dec <- rast("data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_dec.tif")
# 
# # Extract forest type as binary variable (coniferous, deciduous and cloud)
# pixel_training_data_raw$forest_type_cloud <- terra::extract(forest_type_cloud, combined_sample_coords)[,2]
# pixel_training_data_raw$forest_type_con <- terra::extract(forest_type_con, combined_sample_coords)[,2]
# pixel_training_data_raw$forest_type_dec <- terra::extract(forest_type_dec, combined_sample_coords)[,2]
# pixel_training_data_raw$forest_type_cloud[is.nan(pixel_training_data_raw$forest_type_cloud)] <- terra::extract(bornholm_forest_type_cloud, combined_sample_coords)[is.nan(pixel_training_data_raw$forest_type_cloud),2]
# pixel_training_data_raw$forest_type_con[is.nan(pixel_training_data_raw$forest_type_con)] <- terra::extract(bornholm_forest_type_con, combined_sample_coords)[is.nan(pixel_training_data_raw$forest_type_con),2]
# pixel_training_data_raw$forest_type_dec[is.nan(pixel_training_data_raw$forest_type_dec)] <- terra::extract(bornholm_forest_type_dec, combined_sample_coords)[is.nan(pixel_training_data_raw$forest_type_dec),2]
# 
# # Check whether extractions were complete (no more NAs)
# sum(is.na(pixel_training_data_raw$forest_type_cloud))
# sum(is.na(pixel_training_data_raw$forest_type_con))
# sum(is.na(pixel_training_data_raw$forest_type_dec))

# ## Plant available water
# 
# # Load raster
# paw_160cm <- rast("data/predictor_data/plant_available_water/paw_160cm.tif")
# 
# # Extract paw
# pixel_training_data_raw$paw_160cm <- terra::extract(paw_160cm, combined_sample_coords)[,2]

# Bjerreskov et al. 2021 Forest type (coniferous / decidious)

# Load rasters
treetype_bjer_dec <- rast("data/predictor_data/treetype/treetype_bjer_dec.tif")
treetype_bjer_con <- rast("data/predictor_data/treetype/treetype_bjer_con.tif")

# Extract treetype
pixel_training_data_raw$treetype_bjer_dec <- terra::extract(treetype_bjer_dec, combined_sample_coords)[,2]
pixel_training_data_raw$treetype_bjer_con <- terra::extract(treetype_bjer_con, combined_sample_coords)[,2]

# Unload rasters
rm(treetype_bjer_dec)
rm(treetype_bjer_con)

# ## Focal variables - Old vars
# 
# # Load rasters
# a_ptv_focal_3x3 <- rast("data/predictor_data/focal_variables/old_vars/a_ptv_focal_3x3.tif")
# canopy_height_focal_3x3 <- rast("data/predictor_data/focal_variables/old_vars/canopy_height_focal_3x3.tif")
# normalized_z_sd_focal_3x3 <- rast("data/predictor_data/focal_variables/old_vars/normalized_z_sd_focal_3x3.tif")
# 
# # Extract data
# pixel_training_data_raw$a_ptv_focal_3x3 <- terra::extract(a_ptv_focal_3x3, combined_sample_coords)[,2]
# pixel_training_data_raw$canopy_height_focal_3x3 <- terra::extract(canopy_height_focal_3x3, combined_sample_coords)[,2]
# pixel_training_data_raw$normalized_z_sd_focal_3x3 <- terra::extract(normalized_z_sd_focal_3x3, combined_sample_coords)[,2]

## Focal variables - new vars

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
  
pixel_training_data_raw <- full_join(pixel_training_data_raw, focal_vars, by = "sample_id")

## Near surface groundwater(summer)
ns_groundwater_summer_utm32_10m <- rast("data/predictor_data/terraennaert_grundvand_10m/ns_groundwater_summer_utm32_10m.tif")
pixel_training_data_raw$ns_groundwater_summer_utm32_10m <- terra::extract(ns_groundwater_summer_utm32_10m, combined_sample_coords)[,2]
rm(ns_groundwater_summer_utm32_10m)

## Terrons (Peng 2020)
terron_point <- rast("data/predictor_data/terron_maps/terron_point.tif")
pixel_training_data_raw$terron_point <- terra::extract(terron_point, combined_sample_coords)[,2]
rm(terron_point)

## Soil variables from Derek / SustainScapes 

Clay_utm32_10m <- rast("data/predictor_data/soil_layers/Clay_utm32_10m.tif")
Sand_utm32_10m <- rast("data/predictor_data/soil_layers/Sand_utm32_10m.tif")
Soc_utm32_10m <- rast("data/predictor_data/soil_layers/Soc_utm32_10m.tif")
pixel_training_data_raw$Clay_utm32_10m <- terra::extract(Clay_utm32_10m, combined_sample_coords)[,2]
pixel_training_data_raw$Sand_utm32_10m <- terra::extract(Sand_utm32_10m, combined_sample_coords)[,2]
pixel_training_data_raw$Soc_utm32_10m <- terra::extract(Soc_utm32_10m, combined_sample_coords)[,2]
rm(list = c("Clay_utm32_10m", "Sand_utm32_10m", "Soc_utm32_10m"))

# Foliage height diversity
foliage_height_diversity <- rast("data/predictor_data/foliage_height_diversity/foliage_height_diversity.tif")
pixel_training_data_raw$foliage_height_diversity <- terra::extract(foliage_height_diversity, combined_sample_coords)[,2]
rm(foliage_height_diversity)

# Save intermediate backup
save(pixel_training_data_raw, file = "data/training_data/pixel_training.Rda")
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
save(pixel_training_data, file = "data/training_data/pixel_training.Rda")
# pixel_training_data <- pixel_training_data_raw
# load("data/training_data/pixel_training.Rda")

## 6) Quality control and easy random forests ---- 

# Quick PCA
pc <- pixel_training_data %>% 
  filter(!is.na(inland_water_mask)) %>%
  filter(!is.na(sea_mask)) %>%
  st_drop_geometry() %>%
  dplyr::select(-(1:4)) %>%
  dplyr::select(-contains("mask")) %>% 
  # apply(2, function(x) sum(is.na(x)))
  na.omit() %>%
  prcomp(., center = T, scale = T)
forest_value <- pixel_training_data %>% 
  filter(!is.na(inland_water_mask)) %>%
  filter(!is.na(sea_mask)) %>%
  st_drop_geometry() %>%
  dplyr::select(-(2:4)) %>%
  dplyr::select(-contains("mask")) %>% 
  # apply(2, function(x) sum(is.na(x)))
  na.omit() %>% pull(forest_value)
summary(pc)
plot(pc)

fviz_pca_ind(pc, geom.ind = "point", pointshape = 21,
             fill.ind = forest_value,
             palette = "jco", 
             addEllipses = TRUE)

# SPlit data
training_final <- pixel_training_data %>% 
  dplyr::select(-contains("mask")) %>%
  dplyr::select(-biowide_region, -dereks_stratification) %>%
  st_drop_geometry() %>% 
  na.omit() %>%
  mutate(forest_value = factor(forest_value))
colnames(training_final) <- gsub("\\-", "\\.", colnames(training_final))
set.seed(1)
sample_rows <- sample(1:nrow(training_final),
                      round(nrow(training_final) * 0.3), replace = F)
data_valid <- training_final[sample_rows,]
data_train <- training_final[-sample_rows,]


# Try a random forest model
model_formula <- paste0(colnames(data_train)[c(-1,-2)], collapse = " + ")
model_formula <- paste0("forest_value ~ ", model_formula, "", collapse = "")
model_formula <- as.formula(model_formula)
rf_model <- randomForest(model_formula, data = data_train, importance = TRUE)

# predict training dataset
data_train$forest_value_preds <- predict(rf_model, data_train, type ="class")
mean(data_train$forest_value_preds == data_train$forest_value)  
table(data_train$forest_value, data_train$forest_value_preds)  

# predict validaiton dataset
data_valid$forest_value_preds <- predict(rf_model, data_valid, type = "class")

# Validate classification accuracy
mean(data_valid$forest_value_preds == data_valid$forest_value)                    
(conf_matrix <- table(data_valid$forest_value, data_valid$forest_value_preds))
(sens <- conf_matrix[1,1] / sum(conf_matrix[1,]))
(spec <- conf_matrix[2,2] / sum(conf_matrix[2,]))
(tss <- sens + spec - 1)

# check importance of variables
arrange(as.data.frame(importance(rf_model)), desc(MeanDecreaseGini))[1:20,3:4]
