# Train random forests to predict the forest class
# Based on BIODWIDE stratification
# Jakob Assmann j.assmann@bio.au.dk 14 February 2021

# Dependencies
library(caret)
library(tidyverse)
library(sf)
library(doParallel)
library(ranger)

# Load data
load("data/training_data/pixel_training_derek.Rda")
load("data/training_data/pixel_valid_derek.Rda")

# Set pseudo random generator seed
set.seed(24231)

# Prep data frames 
train_data <- pixel_training_derek %>% 
  ungroup() %>%
  dplyr::select(-sample_id, 
                -polygon_source,  
                -biowide_region, 
                -dereks_stratification,
                -cell)
test_data <- pixel_valid_derek %>% 
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
tuneGrid <- expand.grid(mtry = c(2:5),
                        splitrule = c("gini", "extratrees"),
                        min.node.size = c(1, 3, 5, 8)) # Not stumps, range usually between 1-8
rf_fit <- train(forest_value ~ .,
                data = train_data,
                method = "ranger",
                trControl = trainControl(method = "repeatedcv", 
                                         repeats = 5, # Increase later
                                         classProbs = TRUE, 
                                         summaryFunction = twoClassSummary),
                tuneGrid = tuneGrid,
                importance = "permutation",
                metric = "ROC")
rf_fit

# Write out model
save(rf_fit, file = "data/models/final_ranger_model_pixel_derek.Rda")

# Stop cluster
stopCluster(cl)

# Variable importance
summary(rf_fit)
varImp(rf_fit)$importance %>% arrange(desc(Overall))

# Validate on test set
test_preds <- predict(rf_fit, newdata = test_data)
confusionMatrix(data = test_preds, test_data$forest_value)
