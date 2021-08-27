# DK Forest LiDAR - Data preparation script to generate rasters for the vector
# response data
# Jakob Assmann j.assmann@bio.au.dk

# Dependencies
library(terra)
library(raster)
library(parallel)
library(sf)
library(tidyverse)
library(purrr)
library(dplyr)
library(readxl)

# Set seed for pseudo random numbers
set.seed(2308)

# Set path to EcoDes-DK15 data dtm_10m (this will form the base for the rasters)
dtm_10m_folder <- "D:\\Jakob\\dk_nationwide_lidar\\data\\outputs\\dtm_10m\\"

# List dtm_10m files (shell is faster than list.files here)
dtm_10m_files <- shell(paste0("dir /b ", dtm_10m_folder, "*.tif"), intern = T)

# Load geometries for high quality forests
high_quality <- list.files("data/response_variables/high_quality_forests/", 
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
plantations <- read_sf("data/response_variables/low_quality_forests/NST_plantations/LitraPolygoner_region/LitraPolygoner_region.shp")
plantations_meta <- read_excel("data/response_variables/low_quality_forests/NST_plantations/NST  2019 08012019 ber 16012020 til bios_au.xlsx") 
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
low_quality <- read_sf("data/response_variables/low_quality_forests/ikke_p25/ikkeP25_skov.shp") %>%
  list(., plantations) %>%
  map(function(x) select(x, !everything())) %>%
  bind_rows() 

# Save geometries
save(high_quality, file = "data/high_quality_polys.Rda")
save(low_quality, file = "data/low_quality_polys.Rda")

# Load list of Ecodes-DK variabels
ecodes_vrt <- read.table("D:/Jakob/dk_nationwide_lidar/data/outputs/list_of_vrts.txt",
                         stringsAsFactors = F)[,1]
# Remove uneeded variables
ecodes_vrt <- ecodes_vrt[!grepl("point_count",ecodes_vrt)]
ecodes_vrt <- ecodes_vrt[!grepl("point_source",ecodes_vrt)]
ecodes_vrt <- ecodes_vrt[!grepl("building",ecodes_vrt)]

# add full file name
ecodes_vrt <- paste0("D:/Jakob/dk_nationwide_lidar/data/outputs/", ecodes_vrt)

# get sample coordinates
high_quality_sample_coords <- st_sample(high_quality, 100000) %>% 
  st_sf() %>%
  mutate(forest_class = "high") %>% 
  mutate(sample_id = paste0(forest_class, "_", 1:n()))
low_quality_sample_coords <- st_sample(low_quality, 100000) %>% 
  st_sf() %>%
  mutate(forest_class = "low") %>% 
  mutate(sample_id = paste0(forest_class, "_", 1:n()))
save(high_quality_sample_coords, file = "data/high_quality_pixel_sample.Rda")
save(low_quality_sample_coords, file = "data/low_quality_pixel_sample.Rda")
# load("data/high_quality_pixel_sample.Rda")
# load("data/low_quality_pixel_sample.Rda")

# transform geometries to raster CRS
target_crs <- st_crs(raster("D:/Jakob/dk_nationwide_lidar/data/outputs/dtm_10m/dtm_10m_6049_684.tif"))
high_quality_sample_coords <- high_quality_sample_coords %>%
  st_transform(target_crs)
low_quality_sample_coords <- low_quality_sample_coords %>%
  st_transform(target_crs)

# bind and convert to vect
combined_sample_coords <- rbind(high_quality_sample_coords,
                                low_quality_sample_coords) %>%
  as_Spatial() %>%
  vect()
save(combined_sample_coords, file = "data/combined_pixel_sample.Rda")
# load("data/combined_pixel_sample.Rda")

# extract sample across all vrts (except point_source_ids)
extract_fun <- function(index){
  vrt_rast <- rast(ecodes_vrt[index])
  cell_values <- terra::extract(vrt_rast, combined_sample_coords)[,2]
  extractions <- data.frame(
    sample_id = combined_sample_coords$sample_id,
    forest_class = combined_sample_coords$forest_class,
    cell_values = cell_values)
  names(extractions)[3] <- gsub(".*/(.*).vrt", "\\1", ecodes_vrt[index])
  return(extractions)
}

system.time(combined_sample <- lapply(#cl, 
                         seq_along(ecodes_vrt),
                         extract_fun))
save(combined_sample, file = "data/pixel_sample.Rda")
# load("data/pixel_sample.Rda")

# Combine dataframes
pixel_training_data_raw <- combined_sample %>% map(function(x) dplyr::select(x, -forest_class)) %>%
  reduce(full_join, by = "sample_id") %>% full_join(rbind(high_quality_sample_coords, 
                                                          low_quality_sample_coords), 
                                                    ., by = "sample_id")

# Extract forest type data (coniferous vs. broadleaf)
forest_type_cloud <- rast("data/conif_vs_broadleaf/forest_type_cloud.tif")
forest_type_con <- rast("data/conif_vs_broadleaf/forest_type_con.tif")
forest_type_dec <- rast("data/conif_vs_broadleaf/forest_type_dec.tif")
bornholm_forest_type_con <- rast("data/conif_vs_broadleaf/bornholm_forest_type_con.tif")
bornholm_forest_type_dec <- rast("data/conif_vs_broadleaf/bornholm_forest_type_dec.tif")
pixel_training_data_raw$forest_type_cloud <- terra::extract(forest_type_cloud, combined_sample_coords)[,2]
pixel_training_data_raw$forest_type_con <- terra::extract(forest_type_con, combined_sample_coords)[,2]
pixel_training_data_raw$forest_type_dec <- terra::extract(forest_type_dec, combined_sample_coords)[,2]
pixel_training_data_raw$forest_type_cloud[is.na(pixel_training_data_raw$forest_type_cloud)] <- 0
pixel_training_data_raw$forest_type_con[is.na(pixel_training_data_raw$forest_type_con)] <- terra::extract(bornholm_forest_type_con, combined_sample_coords)[is.na(pixel_training_data_raw$forest_type_con),2]
pixel_training_data_raw$forest_type_dec[is.na(pixel_training_data_raw$forest_type_dec)] <- terra::extract(bornholm_forest_type_dec, combined_sample_coords)[is.na(pixel_training_data_raw$forest_type_dec),2]
sum(is.na(pixel_training_data_raw$forest_type_cloud))
sum(is.na(pixel_training_data_raw$forest_type_con))
sum(is.na(pixel_training_data_raw$forest_type_dec))

# Add plant available water
paw_160cm <- rast("data/plant_available_water/paw_160cm.tif")
pixel_training_data_raw$paw_160cm <- terra::extract(paw_160cm, combined_sample_coords)[,2]


# Add focal variables
paw_160cm_focal <- rast("data/focal_variables/a_ptv_focal_3x3.tif")
pixel_training_data_raw$paw_160cm_focal_3x3 <- terra::extract(paw_160cm_focal, combined_sample_coords)[,2]
canopy_height_focal <- rast("data/focal_variables/canopy_height_focal_3x3.tif")
pixel_training_data_raw$canopy_height_focal_3x3 <- terra::extract(canopy_height_focal, combined_sample_coords)[,2]
normalized_zd_sd_focal <- rast("data/focal_variables/normalized_z_sd_focal_3x3.tif")
pixel_training_data_raw$normalized_zd_sd_focal_3x3 <- terra::extract(normalized_zd_sd_focal, combined_sample_coords)[,2]

# Export final data
save(pixel_training_data_raw, file = "data/pixel_sample_combined.Rda")
# load("data/pixel_sample_combined.Rda")

# Add 

# Quick PCA
pc <- pixel_training_data_raw %>% 
  filter(!is.na(inland_water_mask)) %>%
  filter(!is.na(sea_mask)) %>%
  st_drop_geometry() %>%
  dplyr::select(-1, -2) %>%
  dplyr::select(-contains("mask")) %>% 
  # apply(2, function(x) sum(is.na(x)))
  na.omit() %>%
  prcomp(., center = T, scale = T)

summary(pc)
plot(pc)
library(factoextra)
fviz_pca_ind(pc, geom.ind = "point", pointshape = 21,
             fill.ind = pixel_training_data_raw %>% 
               filter(!is.na(inland_water_mask)) %>%
               filter(!is.na(sea_mask)) %>%
               st_drop_geometry() %>% 
               # apply(2, function(x) sum(is.na(x)))
               na.omit() %>% 
               pull(forest_class),
             palette = "jco", 
             addEllipses = TRUE)

library(randomForest)

# SPlit data
training_final <- pixel_training_data_raw %>% 
  dplyr::select(-contains("mask")) %>%
  st_drop_geometry() %>% na.omit() %>%
  mutate(forest_class = factor(forest_class))
colnames(training_final) <- gsub("\\-", "\\.", colnames(training_final))
set.seed(1)
sample_rows <- sample(1:nrow(training_final),
                      round(nrow(training_final) * 0.3), replace = F)
data_valid <- training_final[sample_rows,]
data_train <- training_final[-sample_rows,]


# Try a random forest model
model_formula <- paste0(colnames(data_train)[c(-1,-2)], collapse = " + ")
model_formula <- paste0("forest_class ~ ", model_formula, "", collapse = "")
model_formula <- as.formula(model_formula)
rf_model <- randomForest(model_formula, data = data_train, importance = TRUE)

# predict training dataset
data_train$forest_class_preds <- predict(rf_model, data_train, type ="class")
mean(data_train$forest_class_preds == data_train$forest_class)  
table(data_train$forest_class, data_train$forest_class_preds)  

# predict validaiton dataset
data_valid$forest_class_preds <- predict(rf_model, data_valid, type = "class")

# Validate classification accuracy
mean(data_valid$forest_class_preds == data_valid$forest_class)                    
(conf_matrix <- table(data_valid$forest_class, data_valid$forest_class_preds))
(sens <- conf_matrix[1,1] / sum(conf_matrix[1,]))
(spec <- conf_matrix[2,2] / sum(conf_matrix[2,]))
(tss <- sens + spec - 1)
# check importance of variables
arrange(as.data.frame(importance(rf_model)), desc(MeanDecreaseGini))[1:20,3:4]
