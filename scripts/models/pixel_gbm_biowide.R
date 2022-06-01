# Train gradient boosted regression tree to predict the forest class
# Based on BIODWIDE stratification
# Jakob Assmann j.assmann@bio.au.dk 8 July 2021

# Dependencies
library(caret)
library(tidyverse)
library(sf)
library(doParallel)
library(gbm)

# Load data
load("data/training_data/pixel_training_biowide.Rda")
load("data/training_data/pixel_valid_biowide.Rda")

# Set pseudo random generator seed
set.seed(24231)

# Prep data
train_data <- pixel_training_biowide %>% 
  ungroup() %>%
  dplyr::select(-sample_id, 
                -polygon_source,  
                -biowide_region, 
                -dereks_stratification,
                -cell)
test_data <- pixel_valid_biowide %>% 
  ungroup() %>%
  dplyr::select(-sample_id, 
                -polygon_source,  
                -biowide_region, 
                -dereks_stratification,
                -cell)

# Register parallel cluster
cl <- makePSOCKcluster(46)
#cl <- makePSOCKcluster(16)
registerDoParallel(cl)

# Optimise hyperparameters for boosted regression tree
# 1) Determine optimum number of trees, fixing other parameters
# Note: I do this outside caret as I found caret's interface was not
# very helfpul for efficiently finding the right number of trees.
tuneGrid <- expand.grid(n.trees = seq(300, 10000, 300), # Check range
                        shrinkage = 0.1, # Slow learning rate to start (range 0.001-0.3)
                        n.minobsinnode = 10, # Default value (range 5-15 common)
                        interaction.depth = 3 # Not stumps, range usually between 1-8
       )
gbm_fit <- train(forest_value ~ .,
                 data = train_data,
                 method = "gbm",
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, # Increase later
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary,
                                          selectionFunction = function(
                                            x,
                                            metric, 
                                            maximize){
                                            tolerance(x, metric, tol = 0.5, maximize)
                                          }),
                 tuneGrid = tuneGrid,
                 metric = "ROC")

# Determine optimum number of trees: 
gbm_fit
plot(gbm_fit)
summary(gbm_fit)
# Best model n.trees = 5400

# 2) Tune learning rate
tuneGrid <- expand.grid(n.trees = c(3400, 5400, 8400), # optimal value determined above plus some scope
                        shrinkage = c(0.3, 0.1, 0.05, 0.01, 0.005), 
                        n.minobsinnode = 10, # Default value (range 5-15 common)
                        interaction.depth = 3 # Not stumps, range usually between 1-8
)
gbm_fit <- train(forest_value ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, # Increase later
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary,
                                          selectionFunction = function(
                                            x,
                                            metric, 
                                            maximize){
                                            tolerance(x, metric, tol = 0.5, maximize)
                                          }),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit
plot(gbm_fit)
# n.trees = 5100 and shrinkage = 0.1 seems to be the best effort vs. gain combo.

# 3) Tune tree parameters
tuneGrid <- expand.grid(n.trees = 5400, # optimal value determined above
                        shrinkage = 0.1, # optimal value determined above
                        n.minobsinnode = c(5,10,15), # Default value (range 5-15 common)
                        interaction.depth = c(2,3,5,8) # Not stumps, range usually between 1-8
)
gbm_fit <- train(forest_value ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, # Increase later
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary,
                                          selectionFunction = function(
                                            x,
                                            metric, 
                                            maximize){
                                            tolerance(x, metric, tol = 0.5, maximize)
                                          }),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit # => optimal values (interaction.depth = 8 and n.minobsinnode = 5) at extremes let's try some more
print(gbm_fit)

# Tune tree parameters again
tuneGrid <- expand.grid(n.trees = 5400, # optimal value determined above
                        shrinkage = 0.1, # optimal
                        n.minobsinnode = c(3,5,7), # Default value (range 5-15 common)
                        interaction.depth = c(6,8,10) # Not stumps, range usually between 1-8
)
gbm_fit <- train(forest_value ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, 
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary,
                                          selectionFunction = function(
                                            x,
                                            metric, 
                                            maximize){
                                            tolerance(x, metric, tol = 0.5, maximize)
                                          }),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit 
plot(gbm_fit)
# Final fit n.trees = 5400, interaction.depth = 8, shrinkage = 0.1 and n.minobsinnode = 3.

# Check last final time whether adding more trees with final other parameters
# is worth it. 
tuneGrid <- expand.grid(n.trees = c(5400, 6400, 7400, 8400, 9400), # to be tested
                        shrinkage = 0.1, # optimal
                        n.minobsinnode = 3, # optimal
                        interaction.depth = 8 # optimal
)
gbm_fit <- train(forest_value ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 5, 
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary,
                                          selectionFunction = function(
                                            x,
                                            metric, 
                                            maximize){
                                            tolerance(x, metric, tol = 0.5, maximize)
                                          }),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit 
plot(gbm_fit)
# Looks like the selected parameters are indeed best simplest model.
# n.trees = 5100, interaction.depth = 8, shrinkage = 0.1 and n.minobsinnode = 3.


# Fit final model with 10 fold cross validation
tuneGrid <- expand.grid(n.trees = 5400, # optimal 
                        shrinkage = 0.1, # optimal
                        n.minobsinnode = 3, # optimal
                        interaction.depth = 8 # optimal
)
gbm_fit <- train(forest_value ~ .,
                 data = train_data,
                 method = "gbm",
                 preProc = c("center", "scale"),
                 trControl = trainControl(method = "repeatedcv", 
                                          repeats = 10, 
                                          classProbs = TRUE, 
                                          summaryFunction = twoClassSummary,
                                          verboseIter = TRUE),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
# Save final model
save(gbm_fit, file = "data/models/final_gbm_model_pixel_biowide.Rda")

# Stop cluster
stopCluster(cl)

# Show results
gbm_fit

# Check variable importance
summary(gbm_fit)
sum(summary(gbm_fit)$rel.inf)
# Validate on test set
test_preds <- predict(gbm_fit, newdata = test_data)
confusionMatrix(data = test_preds, test_data$forest_value)

