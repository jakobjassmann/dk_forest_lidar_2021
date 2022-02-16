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

# Rename and subsample for speed
train_data <- pixel_training_derek %>% 
  sample_n(1500) %>%
  #sample_frac(0.5) %>%
  ungroup() %>%
  dplyr::select(-sample_id, -biowide_region, -dereks_stratification) %>%
  dplyr::select(-forest_type_cloud,
                -forest_type_con,
                -heat_load_index,
                -aspect,
                -openness_mean,
                -normalized_z_mean,
                -twi,
                -contains("proportion"),
                -contains("paw"))
test_data <- pixel_valid_derek %>% 
  sample_n(450) %>%
  #sample_frac(0.5) %>%
  ungroup() %>%
  dplyr::select(-sample_id, -biowide_region, -dereks_stratification) %>%
  dplyr::select(-forest_type_cloud,
                -forest_type_con,
                -heat_load_index,
                -aspect,
                -openness_mean,
                -normalized_z_mean,
                -twi,
                -contains("proportion"),
                -contains("paw"))

# Register parallel cluster
cl <- makePSOCKcluster(30)
registerDoParallel(cl)

# Optimise hyperparameters for boosted regression tree
# 1) Determine optimum number of trees, fixing other parameters
tuneGrid <- expand.grid(mtry = c(2:10),
                        splitrule = c("gini", "extratrees"),
                        min.node.size = c(1, 3, 5)) # Not stumps, range usually between 1-8

rf_fit <- train(forest_value ~ .,
                data = train_data,
                method = "ranger",
                preProc = c("center", "scale"),
                trControl = trainControl(method = "repeatedcv", 
                                         repeats = 5, # Increase later
                                         classProbs = TRUE, 
                                         summaryFunction = twoClassSummary),
                tuneGrid = tuneGrid,
                importance = "permutation",
                metric = "ROC")
rf_fit

# Variable importance
summary(rf_fit)
varImp(rf_fit)

save(rf_fit, file = "data/final_ranger_model_pixel_derek.Rda")
