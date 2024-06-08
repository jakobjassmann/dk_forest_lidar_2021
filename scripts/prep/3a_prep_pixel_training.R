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
library(exactextractr)

# Set seed for pseudo random numbers
set.seed(1168)

## 2) Prepare forest training polygons ----

## High quality forests

# Load and prepare geometries
high_quality <- list.files("data/response_data/high_quality_forests/", 
                           ".shp", 
                           recursive = T,
                           full.names = T) %>%
  # Load shapefile, then do first steps of preprocessing
  map(function(shp_file){
    # Read file
    polygons <- read_sf(shp_file)
    # Filter old growth forests if needed (aftaler om natur)
    if(sum("tilskudsor" %in% names(polygons)) > 0){
      polygons <- filter(polygons, tilskudsor == "Privat urørt skov")
    }
    # Remove all auxilliary data
    polygons <- select(polygons, !everything())
    # Clean geometries
    polygons <- polygons %>%
      # Make valid
      st_make_valid() %>%
      # Remove a pixel diagonal from the polygons
      st_buffer(-(sqrt(10^2+10^2))) #%>%
      #st_difference()
    # Add polygons source column based on the file name
    polygons$polygon_source <- gsub(".*/(.*)\\.shp", "\\1", shp_file)
    # Return processed geometries
    return(polygons)
    }) %>%
  # Bind geometries into one sf object
  bind_rows() %>%
  # Clean up polygon source names
  mutate(polygon_source = case_when(
    polygon_source == "skov_kortlaegning_2016_2018" ~ "p15",
    polygon_source == "p25_offentligareal" ~ "p25",
    polygon_source == "aftale_natur_tinglyst" ~ "private_old_growth",
  )) 

## Remove overlap between polygon

# Re arrange order of dataframe: 1) p25, 2) private old growth, 3) p15 forests
high_quality <- high_quality[order(match(high_quality$polygon_source, c("p25", "private_old_growth", "p15"))),]
# add arbitary id
high_quality$id <- paste0(high_quality$polygon_source, "_", 1:nrow(high_quality))

# Check overlap
high_quality$overlaps_with <- high_quality %>%
  st_overlaps()
# Add stats for overlap
high_quality <- high_quality %>%
  # Number of other polygons each polygon overlaps with
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  # Does the polygon overlap with another one at all? T/F
  mutate(overlaps = n_overlaps > 0)
# Print out stats
high_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()
# # A tibble: 12 x 2
# n_overlaps     n
# <int> <int>
#   1          0  8930
# 2          1   301
# 3          2    37
# 4          3    15
# 5          4     7
# 6          5     2
# 7          6     2
# 8          7     1
# 9          8     2
# 10          9     1
# 11         10     1
# 12         13     2

# Clean overlap sequentially in order of sf, starting with the p25 polygons
for(poly_id in high_quality$id[high_quality$overlaps]){
  cat("Cleaning up for", poly_id, "\n")
  # Get polygon geometry
  poly_geo <- high_quality %>% 
    filter(id == poly_id)
  # Retrieve index of overlapping ploygons
  overlaps_with <- high_quality %>% 
    filter(id == poly_id) %>%
    pull(overlaps_with) %>%
    unlist() 
  # Remove from geometry from all other polygons that it overlaps with further 
  # down in the sf objhect
  high_quality[overlaps_with,]$geometry <- overlaps_with %>% 
    slice(high_quality, .) %>%
    split(., 1:nrow(.)) %>%
    map(function(x) st_difference(st_geometry(x), poly_geo)) %>%
    do.call(c, .)
  # Clean up environment
  rm(poly_geo)
  rm(overlaps_with)
}

# Check results (same as above)
high_quality$overlaps_with <- high_quality %>%
  st_overlaps()
high_quality <- high_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
high_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()
# # A tibble: 8 x 2
# n_overlaps     n
# <int> <int>
#   1          0  9116
# 2          1   158
# 3          2    16
# 4          3     6
# 5          4     2
# 6          5     1
# 7          6     1
# 8          7     1

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
# # A tibble: 2 x 2
# n_overlaps     n
# <int> <int>
#   1          0  9299
# 2          1     2
# Good only one overlapping pair remaining. 
test <- high_quality[high_quality$overlaps,]
plot(st_geometry(test[1,]))
plot(st_geometry(test[2,]), add = T, col = "red")
rm(test)
# The private old growth fully contains the p15 forest. Remove the latter
high_quality <- filter(high_quality,
                       !(overlaps & polygon_source == "p15"))

# One last check:
high_quality$overlaps_with <- high_quality %>%
  st_overlaps()
high_quality <- high_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
high_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()

# Check if any empty geometries ar present
sum(st_is_empty(high_quality))

# remove those
high_quality <- filter(high_quality, !st_is_empty(high_quality))

# Check whether any thing else but polygons are present
unique(st_geometry_type(high_quality))

# Great! All cleaned :)

# Remove excess colums
high_quality <- high_quality %>%
  select(polygon_source, geometry)

## Load and prep geometries for low quality forests

# Load ikke p25 geometries
ikke_p25 <- read_sf("data/response_data/low_quality_forests/ikke_p25/ikkeP25_skov.shp") 
# Make geomteries valid
ikke_p25 <- ikke_p25 %>% st_make_valid() 
# Add a polygon source
ikke_p25$polygon_source <- "ikke_p25"
# Remove all unwanted columns
ikke_p25 <- select(ikke_p25, geometry, polygon_source)
# Buffer geometries by a pixel diagonal
ikke_p25 <- st_buffer(ikke_p25, -(sqrt(10^2+10^2)))

# Load plantation geometries
plantations <- read_sf("data/response_data/low_quality_forests/NST_plantations/LitraPolygoner_region/LitraPolygoner_region.shp")
# Load plantation meta data
plantations_meta <- read_excel("data/response_data/low_quality_forests/NST_plantations/NST  2019 08012019 ber 16012020 til bios_au.xlsx") 

# Helper function to classify plantation age bins
# Note, these are given as decadal mid-points, e.g., 5, 15, 25, 35 etc.
sort_into_age_bins <- function(Aldersklasse){
  case_when(Aldersklasse > 0 & Aldersklasse <= 5 ~ "0_to_10",  # Filters decadal mid-point 5
            Aldersklasse > 5 & Aldersklasse <= 25 ~ "10_to_30", # Filters decadal mid-points 15 and 25
            Aldersklasse > 25 & Aldersklasse <= 45 ~ "30_to_50", # Filters decadal mid-points 35 and 45
            Aldersklasse > 45 & Aldersklasse <= 65 ~ "50_to_70", # Filters decadal mid-poinst 55 and 65
            Aldersklasse > 65 & Aldersklasse <= 95 ~ "70_to_100", # Filters decadal mid-points 75, 85 and 95
            TRUE ~ "NA") # Everything else is set to NA
}

# Prepare and clean plantation meta data
plantations_meta <- plantations_meta %>%
  # Remove unwatned columns
  select(Ident, Aldersklasse, `ANV 3`, `Allerede urørt`, Status) %>%
  filter(`ANV 3` == 0) %>% # Keep only ANV 3 values of “0” (remove everything else)
  filter(`Allerede urørt` != "Urørt") %>% # Throw out all rows with label “Urørt"
  filter(Status == "G") %>% # Keep only current plantations (status = “G”) 
  # Remove NA values
  na.omit() %>%
  # Sort into age bins using helper function
  mutate(age_bin = sort_into_age_bins(Aldersklasse)) %>%
  # Remove NAs once more
  filter(age_bin != "NA") %>% # Remove all age classes not in the above categories
  # Group tibble
  group_by(age_bin) %>%
  # Sample each age bin 1000 times
  sample_n(1000)

# Make geometries valid and filter using the meta data
plantations <- plantations %>%
  st_make_valid() %>%
  filter(UNIKID %in% plantations_meta$Ident)
# Set source column
plantations$polygon_source <- "NST_plantations"
# Remove all unwanted columns
plantations <- select(plantations, geometry, polygon_source)
# Buffer geometries by a pixel diagonal
plantations <- st_buffer(plantations, -(sqrt(10^2+10^2)))


# Note: Next we remove overlap between high quality and low quality polygons,
# as well as overlap between low quality polygons

# First we append the low quality polygons to the high quality data
# This is so we can clean from top to bottom again
low_quality <- bind_rows(high_quality,
                         ikke_p25,
                         plantations)

# add a unique polygon id
low_quality$id <- paste0(low_quality$polygon_source, "_", 1:nrow(low_quality))

# Calculate overlap statistics (same as for high quality polygons above)
low_quality$overlaps_with <- low_quality %>%
  st_overlaps()
low_quality <- low_quality %>%
  mutate(n_overlaps = sapply(overlaps_with, length)) %>%
  mutate(overlaps = n_overlaps > 0)
# Check out stats
low_quality %>% st_drop_geometry() %>% group_by(n_overlaps) %>% tally()
# # A tibble: 9 x 2
# n_overlaps     n
# <int> <int>
#   1          0 17651
# 2          1  1551
# 3          2   219
# 4          3    51
# 5          4    15
# 6          5     6
# 7          6     2
# 8          7     3
# 9         11     1

# Remove overlap starting from top to bottom working only on low quality polys
for(poly_id in low_quality$id[low_quality$overlaps]){
  cat("Cleaning up for", poly_id, "\n")
  # Get polygon geometry
  poly_geo <- low_quality %>% 
    filter(id == poly_id)
  # Get index values for polygons that overlap with the target polygon
  overlaps_with <- low_quality %>% 
    filter(id == poly_id) %>%
    pull(overlaps_with) %>%
    unlist() 
  # Remove all polygons high_quality polygons to avoid double cleaning of those
  overlaps_with <- overlaps_with[!(low_quality$polygon_source[overlaps_with] %in% c("p25", "private_old_growth", "p15"))]
  # Remove from geometry from all other polygons that it overlaps with
  if(length(overlaps_with) > 0) {
    low_quality[overlaps_with,]$geometry <- overlaps_with %>% 
    slice(low_quality, .) %>%
    split(., 1:nrow(.)) %>%
    map(function(x){
      difference <- st_difference(st_geometry(x), poly_geo)
      # sf can't handle merging certain types of empty geometries when replacing 
      # the data, st_geometry(st_polygon()) did not work either. 
      # However, creating an empty geometry by over-buffering seemed to work. 
      if(length(difference) == 0) difference <- st_geometry(st_buffer(poly_geo, -10^6))
      return(difference)
    }) %>%
    # Concatenate geometries and assign
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
# A tibble: 8 x 2
# n_overlaps     n
# <int> <int>
#   1          0 17277
# 2          1  1302
# 3          2   132
# 4          3    28
# 5          4     4
# 6          5     3
# 7          6     2
# 8         11     1
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
#   1          0 18747
# 2          1     2
# Only two more forests!
# Let's check them out
test <- low_quality[low_quality$overlaps,]
test
plot(st_geometry(test[1,]))
plot(st_geometry(test[2,]), add = T, col = "red")
rm(test)

# The forests are identical, remove the ikke_p25 forest and finalise the
# low_quality set by removing the high quality forests
low_quality <- filter(low_quality, !(id %in% (low_quality[low_quality$overlaps,] %>% 
                               filter(polygon_source == "ikke_p25") %>%
                               pull(id)))) %>%
  filter(polygon_source %in% c("ikke_p25", "NST_plantations"))
# Remove surplus columns
low_quality <- select(low_quality, polygon_source, geometry)

# Remove empty polygons
low_quality <- low_quality[!st_is_empty(low_quality),]

# Some geometries are MULTILINESTRINGS in the low quality pologyons
unique(st_geometry_type(low_quality))
sum(st_geometry_type(low_quality) == "MULTILINESTRING")
sum(st_geometry_type(low_quality) == "GEOMETRYCOLLECTION")
# 3

# ... remove those MULTILINESTRINGs and "GEOMETRYCOLLECTION"
low_quality <- low_quality[st_geometry_type(low_quality) != "MULTILINESTRING", ]
low_quality <- low_quality[st_geometry_type(low_quality) != "GEOMETRYCOLLECTION", ]

# Final check to make sure there is no overlap
st_overlaps(high_quality) %>% sapply(length) %>% sum()
st_overlaps(low_quality) %>% sapply(length) %>% sum()
st_overlaps(high_quality, low_quality) %>% sapply(length) %>% sum()
# Great! :)

# Save / load geometries
save(high_quality, file = "data/training_data/polygon_geometries/high_quality_polys.Rda")
save(low_quality, file = "data/training_data/polygon_geometries/low_quality_polys.Rda")
# load("data/training_data/polygon_geometries/high_quality_polys.Rda")
# load("data/training_data/polygon_geometries/low_quality_polys.Rda")

## 3) Generate pixel samples from forest polygons ----

# Load EcoDes raster to determine pixels overlapping with polygons
dtm_10m <- rast("F:/JakobAssmann/EcoDes-DK_v1.1.0/dtm_10m/dtm_10m.vrt")

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
                                30000 - nrow(high_quality_sample)) %>%
  bind_rows(high_quality_sample, .)
low_quality_sample <- sample_n(low_quality_pixels,
                                30000 - nrow(low_quality_sample)) %>%
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
save(high_quality_sample, file = "data/training_data/pixel_geometries/high_quality_pixel_sample.Rda")
save(low_quality_sample, file = "data/training_data/pixel_geometries/low_quality_pixel_sample.Rda")
# load("data/training_data/pixel_geometries/high_quality_pixel_sample.Rda")
# load("data/training_data/pixel_geometries/low_quality_pixel_sample.Rda")

# bind into one sf object 
combined_sample_coords <- rbind(high_quality_sample,
                                low_quality_sample)

# Save / load
save(combined_sample_coords, file = "data/training_data/pixel_geometries/combined_pixel_sample.Rda")
# load("data/training_data/pixel_geometries/combined_pixel_sample.Rda")

# Write as shp (not required for parallel processing with terra anymore)
# write_sf(combined_sample_coords, 
#          dsn = "data/training_data/pixel_geometries/combined_pixel_sample.shp")

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
cl <- makeCluster(14)
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

# ## Terrons (Peng 2020)
# 
# # Load raster
# terron_point <- rast("data/predictor_data/terron_maps/terron_point.tif")
# 
# # Extract values
# pixel_training_data_raw$terron_point <- terra::extract(terron_point, combined_sample_coords)[,2]
# 
# # Remove raster
# rm(terron_point)

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
