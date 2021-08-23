# First attempt to train a gradient boosted regression tree to predict the forest class
# for the POLYGONS
# Jakob Assmann j.assmann@bio.au.dk 13 July 2021

# Dependencies
library(raster)
library(tidyverse)
library(sf)
library(caret)
library(gbm)
library(doParallel)

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

# Combined
combined_polys <- rbind(high_quality, 
                        low_quality)

# Transform polygons to raster data CRS
combined_polys <- st_transform(combined_polys,
                               st_crs(raster("D:/Jakob/dk_nationwide_lidar/data/outputs/dtm_10m/dtm_10m_6049_684.tif")))

# Load Training Data
load("data/poly_training.Rda")

# Combine dataframes
poly_training_data <- poly_training_data %>% map(st_drop_geometry) %>% map(function(x) dplyr::select(x, -forest_class)) %>%
  reduce(full_join, by = "forest_id") %>% right_join(combined_polys, ., by = "forest_id")

# Remove - from variable names
names(poly_training_data) <- gsub("-", ".", names(poly_training_data))

# Save tidy poly training data
save(poly_training_data, file = "data/poly_training_tidy.Rda")

# Get rid of NAs
poly_training_data <- na.omit(poly_training_data) %>%
  dplyr::select(-contains("mask"))

# Slit dataset into training and control
set.seed(32423)
index_training <- createDataPartition(
  y = poly_training_data$forest_class,
  p = 0.8,
  list = F
)
train_data <- st_drop_geometry(poly_training_data)[index_training, -1]
test_data <- st_drop_geometry(poly_training_data)[-index_training, -1]

# Register parallel cluster
cl <- makePSOCKcluster(26)
registerDoParallel(cl)

# Optimise hyperparameters for boosted regression tree
# 1) Determine optimum number of trees, fixing other parameters
tuneGrid <- expand.grid(n.trees = seq(200,10000, 300), # Check range
                        shrinkage = 0.1, # Slow learning rate to start (range 0.001-0.3)
                        n.minobsinnode = 10, # Default value (range 5-15 common)
                        interaction.depth = 3 # Not stumps, range usually between 1-8
)
gbm_fit <- train(forest_class ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, # Increase later
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit
# Optimal n.trees = 1100
# 2) Tune learning rate
tuneGrid <- expand.grid(n.trees = 1100, # 1100 optimal value determined above
                        shrinkage = c(0.3, 0.1, 0.05, 0.01, 0.005), 
                        n.minobsinnode = 10, # Default value (range 5-15 common)
                        interaction.depth = 3 # Not stumps, range usually between 1-8
)
gbm_fit <- train(forest_class ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, # Increase later
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit # 0.1 seems to be a good learning rate

# Dig a bit deeper optimising the learning rate
tuneGrid <- expand.grid(n.trees = 1100, # 1100 optimal value determined above
                        shrinkage = seq(0.05, 0.2, 0.025), # search around 0.1
                        n.minobsinnode = 10, # Default value (range 5-15 common)
                        interaction.depth = 3 # Not stumps, range usually between 1-8
)
gbm_fit <- train(forest_class ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, # Increase later
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit # best shrinkage = 0.75

# Tune tree parameters
tuneGrid <- expand.grid(n.trees = 1100, # 1100 optimal value determined above
                        shrinkage = 0.75, # 0.75 optimal
                        n.minobsinnode = c(5,10,15), # Default value (range 5-15 common)
                        interaction.depth = c(2,3,5,8) # Not stumps, range usually between 1-8
)
gbm_fit <- train(forest_class ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, # Increase later
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit # => No improvements let's stick to the previous values


# Fit final model with 10 fold cross validaiton
tuneGrid <- tuneGrid <- expand.grid(n.trees = 1100, 
                                    shrinkage = 0.75,
                                    n.minobsinnode = 10, 
                                    interaction.depth = 3 )

gbm_fit <- train(forest_class ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 10, 
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
# Save final model
save(gbm_fit, file = "data/final_gbm_model_poly.Rda")

# Stop cluster
stopCluster(cl)

# Show results
gbm_fit

# Check variable importance
varImp(gbm_fit)

# Validate on test set
test_preds <- predict(gbm_fit, newdata = test_data)
confusionMatrix(data = test_preds, as.factor(test_data$forest_class))
