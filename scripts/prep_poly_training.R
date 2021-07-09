# DK Forest Polygon Prep
# Jakob Assmann j.assmann@bio.au.dk 7 July 2021

# Dependencies
library(sf)
library(raster)
library(tidyverse)
library(dplyr)
library(purrr)
library(exactextractr)
library(parallel)
library(factoextra)
library(randomForest)

# Load shapefiles
high_quality <- list.files("data/response/high_quality_forests/", 
                           ".shp", 
                           recursive = T,
                           full.names = T) %>%
  map(function(x) {
    sf_object <- read_sf(x)
    sf_object <- sf_object %>%
      select(!everything()) %>%
      mutate(forest_id = paste0(gsub(".*/(.*)\\.shp", "\\1", x),
                                "_", 1:n()),
             forest_class = "high")
    }) %>%
  bind_rows() 
low_quality <- list.files("data/response/low_quality_forests/", 
                          "shp$", 
                          recursive = T,
                          full.names = T) %>%
  map(function(x) {
    sf_object <- read_sf(x)
    sf_object <- sf_object %>%
      select(!everything()) %>%
      mutate(forest_id = paste0(gsub(".*/(.*)\\.shp", "\\1", x),
                                "_", 1:n()),
             forest_class = "low")
  }) %>%
  bind_rows() 

# Create balanced training dataset
combined_polys <- rbind(high_quality, 
                        low_quality[sample(1:nrow(low_quality), 10000),])

# Load list of Ecodes-DK variabels
ecodes_vrt <- read.table("D:/Jakob/dk_nationwide_lidar/data/outputs/list_of_vrts.txt",
                         stringsAsFactors = F)[,1]
# Remove unneded variables
ecodes_vrt <- ecodes_vrt[!grepl("point_count",ecodes_vrt)]
ecodes_vrt <- ecodes_vrt[!grepl("point_source",ecodes_vrt)]
ecodes_vrt <- ecodes_vrt[!grepl("building",ecodes_vrt)]

# add full file name
ecodes_vrt <- paste0("D:/Jakob/dk_nationwide_lidar/data/outputs/", ecodes_vrt)

# Transofrm polygons to raster data CRS
combined_polys <- st_transform(combined_polys,
                               st_crs(raster("D:/Jakob/dk_nationwide_lidar/data/outputs/dtm_10m/dtm_10m_6049_684.tif")))

# Extract mean statistics for polygons using exactextract
cl <- makeCluster(42)
clusterEvalQ(cl, library(exactextractr))
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(sf))
system.time(poly_training_data <- parLapply(cl, 
  ecodes_vrt, 
  function(vrt_file, sample_polys){
    ecodes_raster <- raster(vrt_file)
    sample_polys$mean <- exact_extract(ecodes_raster, sample_polys, fun = "mean")
    sample_polys$sd <- exact_extract(ecodes_raster, sample_polys, fun = "stdev")
    names(sample_polys)[4] <- paste0(gsub(".*/(.*).vrt", "\\1", vrt_file), "_mean")
    names(sample_polys)[5] <- paste0(gsub(".*/(.*).vrt", "\\1", vrt_file), "_sd")
    return(sample_polys)
  },
  combined_polys))
save(poly_training_data, file = "data/poly_training.Rda")
stopCluster(cl)

# Combine dataframes
poly_training_data <- poly_training_data %>% map(st_drop_geometry) %>% map(function(x) dplyr::select(x, -forest_class)) %>%
  reduce(full_join, by = "forest_id") %>% full_join(combined_polys, ., by = "forest_id")

# Quick PCA
pc <- poly_training_data %>% 
  dplyr::select(-contains("mask")) %>%
  st_drop_geometry() %>%
  dplyr::select(-1, -2) %>%
  #apply(2, function(x) sum(is.na(x)))
  na.omit() %>%
  prcomp(., center = T, scale = T)
summary(pc)
plot(pc)
fviz_pca_ind(pc, geom.ind = "point", pointshape = 21,
             fill.ind = poly_training_data %>% 
               dplyr::select(-contains("mask")) %>%
               st_drop_geometry() %>% na.omit() %>% 
               mutate(source = gsub("(.*)(_[0-9]*)$", "\\1", forest_id)) %>%
               pull(source),
             palette = "jco", 
             addEllipses = TRUE)

# SPlit data
training_final <- poly_training_data %>% 
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
importance(rf_model)
