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
tuneGrid <- expand.grid(mtry = c(2:5),
                        splitrule = c("gini", "extratrees"),
                        min.node.size = c(1, 3, 5)) # Not stumps, range usually between 1-8

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

# Variable importnace
summary(rf_fit)
varImp(rf_fit)$importance %>% arrange(desc(Overall))
# Save final model
save(rf_fit, file = "data/models/final_ranger_model_pixel_biowide.Rda")

# Stop Cluster
stopCluster(cl)
