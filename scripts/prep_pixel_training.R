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

# Set path to EcoDes-DK15 data dtm_10m (this will form the base for the rasters)
dtm_10m_folder <- "D:\\Jakob\\dk_nationwide_lidar\\data\\outputs\\dtm_10m\\"

# List drm_10m files (shell is faster than list.files here)
dtm_10m_files <- shell(paste0("dir /b ", dtm_10m_folder, "*.tif"), intern = T)

# Load shapefiles
high_quality <- list.files("data/response/high_quality_forests/", 
                           ".shp", 
                           recursive = T,
                           full.names = T) %>%
  map(read_sf) %>%
  map(function(x) select(x, !everything())) %>%
  bind_rows() 
low_quality <- list.files("data/response/low_quality_forests/", 
                         "shp$", 
                         recursive = T,
                         full.names = T) %>%
  map(read_sf) %>%
  map(function(x) select(x, !everything())) %>%
  bind_rows() 

# Load list of Ecodes-DK variabels
ecodes_vrt <- read.table("D:/Jakob/dk_nationwide_lidar/data/outputs/list_of_vrts.txt",
                         stringsAsFactors = F)[,1]
# Remove unneded variables
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
forest_type <- rast("data/conif_vs_broadleaf/dk_forest_con_vs_dec.tif")
pixel_training_data_raw$forest_type <- terra::extract(forest_type, combined_sample_coords)[,2]

# Export final data
save(pixel_training_data_raw, file = "data/pixel_sample_combined.Rda")

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
