# Sample size sensitivity analysis for the random forest models 
# Based on BIODWIDE stratification and best performing hyperparameters for
# n = 30k

# Jakob Assmann jakob.assmann@uzh.ch 1 June 2024

# Dependencies
library(caret)
library(tidyverse)
library(sf)
library(pbapply)
library(ranger)
library(gbm)
library(ggplot2)
library(cowplot)

# Set pseudo random generator seed
set.seed(240601)

# Load data for sensitivity analysis (max n approx. 130k)
load("data/training_data/pixel_training_biowide_sensitivity.Rda")
load("data/training_data/pixel_valid_biowide_sensitivity.Rda")
load("data/training_data/pixel_training_derek_sensitivity.Rda")
load("data/training_data/pixel_valid_derek_sensitivity.Rda")

# Prep data frames biowide data
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

# Helper function for sub-sample training and validation data to a total target
get_sub_train <- function(n_target){
  sample_n(train_data, round(n_target * 0.8))
}
get_sub_test <- function(n_target){
  sample_n(test_data, round(n_target * 0.2))
}
# Quick test
nrow(get_sub_train(100000))
nrow(get_sub_test(100000))

# Set target sampling sizes to test
n_samples <- c(10000, 20000, 30000, 40000, 50000, 75000, 100000, 125000)

# Ranger training function 
train_ranger <- function(n){
  # Train model on subset of data and time
  train_data_sub <- get_sub_train(n)
  time_training <- system.time(
    rf_fit <- train(forest_value ~ .,
                    data = train_data_sub,
                    method = "ranger",
                    trControl = trainControl(method = "none",
                                             classProbs = TRUE, 
                                             summaryFunction = twoClassSummary),
                    tuneGrid = hyper_params,
                    importance = "permutation",
                    metric = "ROC")
  )
  
  # Validate on test set
  test_data_sub <- get_sub_test(n) 
  test_preds <- predict(rf_fit, newdata = test_data_sub)
  cf_mat <- confusionMatrix(data = test_preds, test_data_sub$forest_value)
  
  # Return confusion matrix and time3
  return(list(cf_mat, time_training))
}

# GBM training function 
train_gbm <- function(n){
  # Train model on subset of data and time
  train_data_sub <- get_sub_train(n)
  time_training <- system.time(
    gbm_fit <- train(forest_value ~ .,
                     data = train_data_sub,
                     method = "gbm",
                     preProc = c("center", "scale"),
                     trControl = trainControl(method = "none", 
                                              classProbs = TRUE, 
                                              summaryFunction = twoClassSummary),
                     tuneGrid = hyper_params,
                     metric = "ROC")
  )
  
  # Validate on test set
  test_data_sub <- get_sub_test(n) 
  test_preds <- predict(gbm_fit, newdata = test_data_sub)
  cf_mat <- confusionMatrix(data = test_preds, test_data_sub$forest_value)
  
  # Return confusion matrix and time3
  return(list(cf_mat, time_training))
}

# Helper function to visualize results
vis_results <- function(results_df, model_name){
  ranger_biowide_plot <- plot_grid(
    ggplot(results_df, aes(x = n_samples, y = accuracy)) +
      geom_line() +
      geom_point() +
      annotate("point", x = 30000, y = results_df$accuracy[3],
               colour = "red") +
      annotate("text", x = 10000, y = Inf, vjust = 1.5, hjust = 0,
               label = model_name, size = 14 / .pt) +
      scale_x_continuous(limits = c(0, 125000)) +
      #scale_y_continuous(limits = c(0.75, 0.9)) +
      labs(x = "Sample size", y = "Accuracy") +
      theme_cowplot(),
    ggplot(results_df, aes(x = n_samples, y = speed)) +
      geom_line() +
      geom_point() +
      annotate("point", x = 30000, y = results_df$speed[3],
               colour = "red") +
      annotate("text", x = 10000, y = Inf, vjust = 1.5, hjust = 0,
               label = model_name, size = 14 / .pt) +
      scale_x_continuous(limits = c(0, 125000)) +
      # scale_y_continuous(limits = c(0, 500)) +
      labs(x = "Sample size", y = "Training time (s)") +
      theme_cowplot()
  )
}

## 1) Ranger Biowide

# Set best performing hyper parameters based on n = 30k 
hyper_params <- data.frame(mtry =5,
                           min.node.size =1,
                           splitrule	= "gini")

# Map over sampling intervals
sens_ranger_biowide <- pblapply(n_samples, train_ranger)

# Generate summary data frame
results_ranger_biowide <- data.frame(
  n_samples = n_samples,
  accuracy = sapply(sens_ranger_biowide, function(x) x[[1]]$overall[1]),
  speed = sapply(sens_ranger_biowide, function(x) x[[2]][1])
)

plot_ranger_biowide <- vis_results(results_ranger_biowide,
                                   "Random Forest\nBiowide")

## 2) GBM biowide

# Set best performing hyper parameters based on n = 30k 
hyper_params <- data.frame(n.trees = 5400,
                           shrinkage = 0.1,
                           n.minobsinnode = 3,
                           interaction.depth = 8)

# Map over sampling intervals
sens_gbm_biowide <- pblapply(n_samples, train_gbm)

# Generate summary data frame
results_gbm_biowide <- data.frame(
  n_samples = n_samples,
  accuracy = sapply(sens_gbm_biowide, function(x) x[[1]]$overall[1]),
  speed = sapply(sens_gbm_biowide, function(x) x[[2]][1])
)

plot_gbm_biowide <- vis_results(results_gbm_biowide,
                                   "gbm Biowide")

## 3) Prepare Sustainscapes data

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


## 4) Ranger Sustainscapes

# Set best performing hyper parameters based on n = 30k 
hyper_params <- data.frame(mtry =5,
                           min.node.size =1,
                           splitrule	= "gini")

# Map over sampling intervals
sens_ranger_derek <- pblapply(n_samples, train_ranger)

# Generate summary data frame
results_ranger_derek <- data.frame(
  n_samples = n_samples,
  accuracy = sapply(sens_ranger_derek, function(x) x[[1]]$overall[1]),
  speed = sapply(sens_ranger_derek, function(x) x[[2]][1])
)

plot_ranger_derek <- vis_results(results_ranger_derek,
                                   "Random Forest\nSustainScapes")

## 5) GBM Sustainscapes

# Set best performing hyper parameters based on n = 30k 
hyper_params <- data.frame(n.trees = 5700,
                           shrinkage = 0.1,
                           n.minobsinnode = 3,
                           interaction.depth = 8)

# Map over sampling intervals
sens_gbm_derek <- pblapply(n_samples, train_gbm)

# Generate summary data frame
results_gbm_derek <- data.frame(
  n_samples = n_samples,
  accuracy = sapply(sens_gbm_derek, function(x) x[[1]]$overall[1]),
  speed = sapply(sens_gbm_derek, function(x) x[[2]][1])
)

plot_gbm_derek <- vis_results(results_gbm_derek,
                                "gbm SustainScapes")

## 6) Generate summary plots
(ranger_sensitivity <- plot_grid(plot_ranger_biowide,
                                plot_ranger_derek,
                                nrow = 2,
                                align = "hv"))

(gbm_sensitivity <- plot_grid(plot_gbm_biowide,
                                plot_gbm_derek,
                                nrow = 2,
                                align = "hv"))

# Save results for back up 
save(results_ranger_biowide, 
     results_ranger_derek, 
     results_gbm_derek,
     results_gbm_biowide,
     file = "data/models/sensitivity_analysis.Rda")
