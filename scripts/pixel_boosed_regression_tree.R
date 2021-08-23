# First attempt to train a gradient boosted regression tree to predict the forest class
# Jakob Assmann j.assmann@bio.au.dk 8 July 2021

# Dependencies
library(caret)
library(tidyverse)
library(sf)
library(doParallel)
library(gbm)

# Load data
load("data/pixel_sample_combined.Rda")

# # Reshape data 
# combined_sample <- combined_sample %>% map(function(x) dplyr::select(x, -forest_class)) %>%
#   reduce(full_join, by = "sample_id") %>% full_join(combined_sample[[1]][,1:2], 
#                                                     ., by = "sample_id")

# Filter out masked data and na values   
combined_sample <- pixel_training_data_raw %>% 
  st_drop_geometry() %>%
  filter(!is.na(inland_water_mask)) %>%
  filter(!is.na(sea_mask)) %>%
  select(-contains("mask")) %>%
  mutate(forest_class = factor(forest_class)) %>%
  na.omit()

# Remove - from variable names
names(combined_sample) <- gsub("-", ".", names(combined_sample))

# Slit dataset into training and control
set.seed(32423)
index_training <- createDataPartition(
  y = combined_sample$forest_class,
  p = 0.8,
  list = F
)
train_data <- combined_sample[index_training, -which(names(combined_sample) == "sample_id")]
test_data <- combined_sample[-index_training, -which(names(combined_sample) == "sample_id")]

# Register parallel cluster
cl <- makePSOCKcluster(32)
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
gbm_fit # 0.05 seems to be a good learning rate

# Dig a bit deeper optimising the learning rate
tuneGrid <- expand.grid(n.trees = 1100, # 1100 optimal value determined above
                        shrinkage = seq(0.01, 0.1, 0.01), # search around 0.05
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
gbm_fit # best shrinkage = 0.06

# Tune tree parameters
tuneGrid <- expand.grid(n.trees = 1100, # 1100 optimal value determined above
                        shrinkage = 0.06, # 0.06 optimal
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
gbm_fit # => optimal values (15 and 8) at upper extreme let's try some more

# Tune tree parameters again
tuneGrid <- expand.grid(n.trees = 1100, # 1100 optimal value determined above
                        shrinkage = 0.06, # 0.06 optimal
                        n.minobsinnode = c(10,15,20), # Default value (range 5-15 common)
                        interaction.depth = c(6,8,10) # Not stumps, range usually between 1-8
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
gbm_fit #15 n.minobsinode is optimal, 
# while additional gain from increasing interactoin.depth is minmal leave at 10


# Fit final model with 10 fold cross validaiton
tuneGrid <- expand.grid(n.trees = 1100, # 1100 optimal value determined above
                        shrinkage = 0.06, # 0.06 optimal
                        n.minobsinnode = 15, # Default value (range 5-15 common)
                        interaction.depth = 10 # Not stumps, range usually between 1-8
)
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
save(gbm_fit, file = "data/final_gbm_model_pixel.Rda")

# Stop cluster
stopCluster(cl)

# Show results
gbm_fit

# Check variable importance
summary(gbm_fit)

# Validate on test set
test_preds <- predict(gbm_fit, newdata = test_data)
confusionMatrix(data = test_preds, test_data$forest_class)

