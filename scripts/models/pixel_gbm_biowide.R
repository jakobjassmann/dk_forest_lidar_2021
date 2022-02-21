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

# Rename and subsample for speed
train_data <- pixel_training_biowide %>% 
  sample_n(1500) %>%
  #sample_frac(0.5) %>%
  ungroup() %>%
  dplyr::select(-sample_id, -biowide_region, -dereks_stratification) 
test_data <- pixel_valid_biowide %>% 
  sample_n(450) %>%
  #sample_frac(0.5) %>%
  ungroup() %>%
  dplyr::select(-sample_id, -biowide_region, -dereks_stratification) 

# Register parallel cluster
cl <- makePSOCKcluster(30)
registerDoParallel(cl)

# Optimise hyperparameters for boosted regression tree
# 1) Determine optimum number of trees, fixing other parameters
tuneGrid <- expand.grid(n.trees = seq(300, 15000, 300), # Check range
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
                                         summaryFunction = twoClassSummary),
                 tuneGrid = tuneGrid,
                 metric = "ROC")
gbm_fit
# n.trees  ROC        Sens       Spec     
# 300    0.8211565  0.8072020  0.6750418
# 600    0.8301944  0.8102908  0.6871564
# 900    0.8347270  0.8113901  0.6939534
# 1200    0.8375925  0.8114438  0.6984440
# 1500    0.8396572  0.8106447  0.7022348
# 1800    0.8411953  0.8105321  0.7055166
# 2100    0.8422804  0.8097545  0.7083001
# 2400    0.8433101  0.8091647  0.7101504
# 2700    0.8440904  0.8086337  0.7124302
# 3000    0.8447755  0.8074594  0.7140261
# 3300    0.8452758  0.8080117  0.7153727
# 3600    0.8457087  0.8070518  0.7165868
# ROC stabalises very quickly. Best Sens / Spec at around n.trees = 1500/1800
# 2) Tune learning rate
tuneGrid <- expand.grid(n.trees = 1200, # 1100 optimal value determined above
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
gbm_fit <- train(forest_value ~ .,
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
gbm_fit <- train(forest_value ~ .,
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
tuneGrid <- expand.grid(n.trees = 1500, # 1100 optimal value determined above
                        shrinkage = 0.06, # 0.06 optimal
                        n.minobsinnode = c(10,15,20), # Default value (range 5-15 common)
                        interaction.depth = c(6,8,10) # Not stumps, range usually between 1-8
)
gbm_fit <- train(forest_value ~ .,
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
tuneGrid <- expand.grid(n.trees = 1500, # 1100 optimal value determined above
                        shrinkage = 0.06, # 0.06 optimal
                        n.minobsinnode = 15, # Default value (range 5-15 common)
                        interaction.depth = 10 # Not stumps, range usually between 1-8
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

