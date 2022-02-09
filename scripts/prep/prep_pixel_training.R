# DK Forest LiDAR - Data preparation script to generate training data
# Jakob Assmann j.assmann@bio.au.dk 4 February 2022

# 1) Housekeeping ----

# Dependencies
library(terra)
library(raster)
library(parallel)
library(pblapply)
library(sf)
library(tidyverse)
library(purrr)
library(dplyr)
library(readxl)
library(factoextra)
library(ggplot2)
library(randomforest)

# Set seed for pseudo random numbers
set.seed(2308)

## 2) Prepare forest training polygons ----

# Load geometries for high quality forests
high_quality <- list.files("data/response_data/high_quality_forests/", 
                           ".shp", 
                           recursive = T,
                           full.names = T) %>%
  map(read_sf) %>%
  map(function(x){
    if(sum("tilskudsor" %in% names(x)) > 0) x <- filter(x, tilskudsor == "Privat urørt skov")
    select(x, !everything())
    }) %>%
  bind_rows() 

## Load and prep geometries for low quality forests
# Plantations geometries and meta data
plantations <- read_sf("data/response_data/low_quality_forests/NST_plantations/LitraPolygoner_region/LitraPolygoner_region.shp")
plantations_meta <- read_excel("data/response_data/low_quality_forests/NST_plantations/NST  2019 08012019 ber 16012020 til bios_au.xlsx") 
# Helper function to classify age bins
sort_into_age_bins <- function(Aldersklasse){
  case_when(Aldersklasse > 0 & Aldersklasse <= 10 ~ "0_to_10",
            Aldersklasse > 10 & Aldersklasse <= 25 ~ "10_to_25",
            Aldersklasse > 25 & Aldersklasse <= 50 ~ "25_to_50",
            Aldersklasse > 50 & Aldersklasse <= 75 ~ "50_to_75",
            Aldersklasse > 75 & Aldersklasse <= 100 ~ "75_to_100",
            TRUE ~ "NA") 
}
# Data cleaning
plantations_meta <- plantations_meta %>%
  select(Ident, Aldersklasse, `ANV 4`, `Allerede urørt`, Status) %>%
  filter(`ANV 4` != 1) %>%
  filter(`Allerede urørt` != "Urørt") %>%
  filter(Status != "H") %>%
  na.omit() %>%
  mutate(age_bin = sort_into_age_bins(Aldersklasse)) %>%
  filter(age_bin != "NA") %>% 
  group_by(age_bin) %>%
  sample_n(1000)
# Filter geometries
plantations <- plantations %>%
  filter(UNIKID %in% plantations_meta$Ident)
# Add ikke p25 geometries
low_quality <- read_sf("data/response_data/low_quality_forests/ikke_p25/ikkeP25_skov.shp") %>%
  list(., plantations) %>%
  map(function(x) select(x, !everything())) %>%
  bind_rows() 

# Save / load geometries
save(high_quality, file = "data/training_data/polygon_geometries/high_quality_polys.Rda")
save(low_quality, file = "data/training_data/polygon_geometries/low_quality_polys.Rda")
# load("data/training_data/polygon_geometries/high_quality_polys.Rda")
# load("data/training_data/polygon_geometries/low_quality_polys.Rda")

## 3) Generate pixel samples from forest polygons ----

# get sample coordinates
high_quality_sample_coords <- st_sample(high_quality, 100000) %>% 
  st_sf() %>%
  mutate(forest_value = "high") %>% 
  mutate(sample_id = paste0(forest_value, "_", 1:n()))
low_quality_sample_coords <- st_sample(low_quality, 100000) %>% 
  st_sf() %>%
  mutate(forest_value = "low") %>% 
  mutate(sample_id = paste0(forest_value, "_", 1:n()))

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

# Write as shp (required for parallel processing with terra later)
write_sf(combined_sample_coords, 
         dsn = "data/training_data/pixel_geometries/combined_pixel_sample.shp")

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
    sample_id = combined_sample_coords$sampl_d,
    forest_value = combined_sample_coords$frst_vl,
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
# load("data/ecodes_pixel_sample_temp.Rda")

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

## Forest type (coniferous vs. broadleaf)

# Load rasters
forest_type_cloud <- rast("data/predictor_data/conif_vs_broadleaf/forest_type_cloud.tif")
forest_type_con <- rast("data/predictor_data/conif_vs_broadleaf/forest_type_con.tif")
forest_type_dec <- rast("data/predictor_data/conif_vs_broadleaf/forest_type_dec.tif")
bornholm_forest_type_cloud <- rast("data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_cloud.tif")
bornholm_forest_type_con <- rast("data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_con.tif")
bornholm_forest_type_dec <- rast("data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_dec.tif")

# Extract forest type as binary variable (coniferous, deciduous and cloud)
pixel_training_data_raw$forest_type_cloud <- terra::extract(forest_type_cloud, combined_sample_coords)[,2]
pixel_training_data_raw$forest_type_con <- terra::extract(forest_type_con, combined_sample_coords)[,2]
pixel_training_data_raw$forest_type_dec <- terra::extract(forest_type_dec, combined_sample_coords)[,2]
pixel_training_data_raw$forest_type_cloud[is.nan(pixel_training_data_raw$forest_type_cloud)] <- terra::extract(bornholm_forest_type_cloud, combined_sample_coords)[is.nan(pixel_training_data_raw$forest_type_cloud),2]
pixel_training_data_raw$forest_type_con[is.nan(pixel_training_data_raw$forest_type_con)] <- terra::extract(bornholm_forest_type_con, combined_sample_coords)[is.nan(pixel_training_data_raw$forest_type_con),2]
pixel_training_data_raw$forest_type_dec[is.nan(pixel_training_data_raw$forest_type_dec)] <- terra::extract(bornholm_forest_type_dec, combined_sample_coords)[is.nan(pixel_training_data_raw$forest_type_dec),2]

# Check whether extractions were complete (no more NAs)
sum(is.na(pixel_training_data_raw$forest_type_cloud))
sum(is.na(pixel_training_data_raw$forest_type_con))
sum(is.na(pixel_training_data_raw$forest_type_dec))

## Plant available water

# Load raster
paw_160cm <- rast("data/predictor_data/plant_available_water/paw_160cm.tif")

# Extract paw
pixel_training_data_raw$paw_160cm <- terra::extract(paw_160cm, combined_sample_coords)[,2]

## Aggregated variables

# Load rasters
paw_160cm_focal <- rast("data/predictor_data/focal_variables/a_ptv_focal_3x3.tif")
canopy_height_focal <- rast("data/predictor_data/focal_variables/canopy_height_focal_3x3.tif")
normalized_zd_sd_focal <- rast("data/predictor_data/focal_variables/normalized_z_sd_focal_3x3.tif")

# Extract data
pixel_training_data_raw$paw_160cm_focal_3x3 <- terra::extract(paw_160cm_focal, combined_sample_coords)[,2]
pixel_training_data_raw$canopy_height_focal_3x3 <- terra::extract(canopy_height_focal, combined_sample_coords)[,2]
pixel_training_data_raw$normalized_zd_sd_focal_3x3 <- terra::extract(normalized_zd_sd_focal, combined_sample_coords)[,2]


## Near surface groundwater(summer)

# Load raster
ns_groundwater_summer <- rast("data/predictor_data/terraennaert_grundvand_10m/summer_predict.tif")
# Extract data
pixel_training_data_raw$ns_groundwater_summer <- terra::extract(ns_groundwater_summer, combined_sample_coords)[,2]

## Terrons (Peng 2020)

# Load raster 
terrons <- rast("data/predictor_data/terron_maps/terron_point.tif")
# Extract data
pixel_training_data_raw$terrons <- terra::extract(terrons, combined_sample_coords)[,2]

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
                                biowide_region,
                                dereks_stratification)
save(pixel_training_data, file = "data/training_data/pixel_training.Rda")
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
