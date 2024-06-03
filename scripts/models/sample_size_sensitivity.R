# Sample size sensitivity analysis for the random forest models 
# Based on BIODWIDE stratification and best performing hyperparameters for
# n = 30k

# Jakob Assmann jakob.assmann@uzh.ch 1 June 2024

# Dependencies
library(caret)
library(tidyverse)
library(sf)
library(pbapply)
library(parallel)
library(ranger)
library(gbm)
library(ggplot2)
library(cowplot)
library(gt)

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
                -cell)
test_data <- pixel_valid_biowide %>% 
  ungroup() %>%
  dplyr::select(-sample_id, 
                -polygon_source,
                -cell)

# Helper function for sub-sample training and validation data to a total target
get_sub_train <- function(n_target, strat_col){
  strat_col <- c(strat_col, "forest_value")
  frac <- 0.8 * 2 * n_target / nrow(train_data)
  train_data %>%
    group_by(across(all_of(!!strat_col))) %>% 
  sample_frac(frac) %>%
    ungroup() %>%
    select(-biowide_region, 
           -dereks_stratification)
}
get_sub_test <- function(n_target, strat_col){
  strat_col <- c(strat_col, "forest_value")
  frac <- 0.2 * 2 * n_target / nrow(test_data) 
  test_data %>%
    group_by(across(all_of(!!strat_col))) %>% 
    sample_frac(frac) %>%
    ungroup() %>%
    select(-biowide_region, 
           -dereks_stratification)
}
# Quick test
nrow(get_sub_train(30000, "biowide_region"))
nrow(get_sub_test(30000, "biowide_region"))

# get_sub_train(30000, "biowide_region") %>% group_by(biowide_region) %>% tally() %>% 
#   mutate(frac = n / group_by(train_data, biowide_region) %>% tally() %>% pull(n),
#          target_frac = 0.8 * 30000 / nrow(train_data))
# get_sub_test(30000, "biowide_region") %>% group_by(biowide_region) %>% tally() %>% 
#   mutate(frac = n / group_by(test_data, biowide_region) %>% tally() %>% pull(n),
#          target_frac = 0.2 * 30000 / nrow(test_data))

get_sub_train(10000, "biowide_region") %>% group_by(forest_value) %>% tally()
get_sub_test(10000, "biowide_region") %>% group_by(forest_value) %>% tally()

# Set target sampling sizes to test
n_samples <- c(10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000)

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
      annotate("text", x = 0, y = Inf, vjust = 1.5, hjust = 0,
               label = model_name, size = 14 / .pt) +
      scale_x_continuous(limits = c(0, 125000)) +
      scale_y_continuous(limits = c(0.75, 0.9)) +
      labs(x = "Sample size", y = "Accuracy") +
      theme_cowplot(),
    ggplot(results_df, aes(x = n_samples, y = speed)) +
      geom_line() +
      geom_point() +
      annotate("point", x = 30000, y = results_df$speed[3],
               colour = "red") +
      annotate("text", x = 0, y = Inf, vjust = 1.5, hjust = 0,
               label = model_name, size = 14 / .pt) +
      scale_x_continuous(limits = c(0, 125000)) +
      scale_y_continuous(limits = c(0, 2500)) +
      labs(x = "Sample size", y = "Training time (s)") +
      theme_cowplot()
  )
}

# Set up cluster
cl <- makeCluster(8)
clusterEvalQ(cl, library(ranger))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(caret))
clusterEvalQ(cl, library(ranger))
clusterEvalQ(cl, library(gbm))
clusterExport(cl, varlist = list("train_data", "test_data", 
                                 "get_sub_train", "get_sub_test",
                                 "train_ranger", "train_gbm"))

## 1) Ranger Biowide

# Set best performing hyper parameters based on n = 30k 
hyper_params <- data.frame(mtry =5,
                           min.node.size =1,
                           splitrule	= "gini")
clusterEvalQ(cl, rm(hyper_params))
clusterExport(cl, "hyper_params")

# Map over sampling intervals
sens_ranger_biowide <- pblapply(n_samples, train_ranger, cl = cl)

# Generate summary data frame
results_ranger_biowide <- data.frame(
  n_samples = n_samples,
  accuracy = sapply(sens_ranger_biowide, function(x) x[[1]]$overall[1]),
  speed = sapply(sens_ranger_biowide, function(x) x[[2]][1])
)

plot_ranger_biowide <- vis_results(results_ranger_biowide,
                                   "Random Forest Biowide")

## 2) GBM biowide

# Set best performing hyper parameters based on n = 30k 
hyper_params <- data.frame(n.trees = 5400,
                           shrinkage = 0.1,
                           n.minobsinnode = 3,
                           interaction.depth = 8)
clusterEvalQ(cl, rm(hyper_params))
clusterExport(cl, "hyper_params")

# Map over sampling intervals
sens_gbm_biowide <- pblapply(n_samples, train_gbm, cl = cl)

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
clusterEvalQ(cl, rm(list = c("train_data", "test_data")))
clusterExport(cl, list("train_data", "test_data"))

## 4) Ranger Sustainscapes

# Set best performing hyper parameters based on n = 30k 
hyper_params <- data.frame(mtry =5,
                           min.node.size =1,
                           splitrule	= "gini")
clusterEvalQ(cl, rm(hyper_params))
clusterExport(cl, "hyper_params")

# Map over sampling intervals
sens_ranger_derek <- pblapply(n_samples, train_ranger, cl = cl)

# Generate summary data frame
results_ranger_derek <- data.frame(
  n_samples = n_samples,
  accuracy = sapply(sens_ranger_derek, function(x) x[[1]]$overall[1]),
  speed = sapply(sens_ranger_derek, function(x) x[[2]][1])
)

plot_ranger_derek <- vis_results(results_ranger_derek,
                                   "Random Forest SustainScapes")

## 5) GBM Sustainscapes

# Set best performing hyper parameters based on n = 30k 
hyper_params <- data.frame(n.trees = 5700,
                           shrinkage = 0.1,
                           n.minobsinnode = 3,
                           interaction.depth = 8)
clusterEvalQ(cl, rm(hyper_params))
clusterExport(cl, "hyper_params")

# Map over sampling intervals
sens_gbm_derek <- pblapply(n_samples, train_gbm, cl = cl)

# Generate summary data frame
results_gbm_derek <- data.frame(
  n_samples = n_samples,
  accuracy = sapply(sens_gbm_derek, function(x) x[[1]]$overall[1]),
  speed = sapply(sens_gbm_derek, function(x) x[[2]][1])
)

plot_gbm_derek <- vis_results(results_gbm_derek,
                                "gbm SustainScapes")

stopCluster(cl)

## 6) Generate summary plots
(ranger_sensitivity <- plot_grid(plot_ranger_biowide,
                                 plot_ranger_derek,
                                 nrow = 2,
                                 align = "hv",
                                 labels = "auto"))
save_plot("docs/figures/sample_sensitivity_ranger.png", 
          ranger_sensitivity,
          nrow = 2,
          ncol = 2,
          bg = "white")
(gbm_sensitivity <- plot_grid(plot_gbm_biowide,
                              plot_gbm_derek,
                              nrow = 2,
                              align = "hv",
                              labels = "auto"))
save_plot("docs/figures/sample_sensitivity_gbm.png", 
          gbm_sensitivity,
          nrow = 2,
          ncol = 2,
          bg = "white")

# Export table as html
bind_rows(results_ranger_biowide %>%
            mutate(model_type = "Random Forest",
                   stratification = "Biowide") %>%
            relocate(c(4,5)),
          results_ranger_derek %>%
            mutate(model_type = "Random Forest",
                   stratification = "SustainScapes") %>%
            relocate(c(4,5)),
          results_gbm_biowide %>%
            mutate(model_type = "gbm",
                   stratification = "Biowide") %>%
            relocate(c(4,5)),
          results_gbm_derek %>%
            mutate(model_type = "gbm",
                   stratification = "SustainScapes") %>%
            relocate(c(4,5))) %>%
  mutate(accuracy = round(accuracy, 2),
         speed = round(speed)) %>%
  select(`Model type` = model_type,
         `Stratification` = stratification,
         `Samples (n)` = n_samples,
         `Accuracy` = accuracy,
         `Speed (s)` = speed) %>%
  group_by(`Model type`, `Stratification`) %>%
  group_map(function(x, ...) x %>%
              gt(rowname_col = "Samples (n)") %>%
              as_raw_html()) %>% 
  data.frame() %>%
  setNames(c("Random Forest Biowide",
             "Random Forest SustainScapes",
             "gbm Biowide",
             "gbm SunstainScapes")) %>%
  gt() %>%
  fmt_markdown(columns = TRUE) %>%
  gtsave("docs/training_size_sensitivity.html")

# Save results for back up 
save(results_ranger_biowide, 
     results_ranger_derek, 
     results_gbm_derek,
     results_gbm_biowide,
     file = "data/models/sensitivity_analysis.Rda")

## 7) Assess nearest neighbours
load("data/training_data/pixel_training_sensitivity.Rda")
pixel_training_data <- pixel_training_data %>% 
  select(-contains("mask")) %>%
  select(-contains("250")) %>%
  select(-contains("mean_110")) %>%
  select(-contains("proportion")) %>%
  mutate(forest_value = factor(forest_value)) %>%
  na.omit()
st_crs(pixel_training_data)

# Get a feeling for how close points are to eachother
sample_1000 <- sample_n(pixel_training_data, 1000)
nearest_feat <- st_nearest_feature(sample_1000)
dist_nf <- sapply(1:1000, function(x){
  st_distance(sample_1000[x,], sample_1000[nearest_feat[x],])
})
mean(dist_nf)
max(dist_nf)
min(dist_nf)

# define function to calculate distances to nearest nieghbour based on buffering
dist_nearest_neighbour <- function(current_id, buffer, sf_obj){
  neighbours <-  sf_obj %>%
    filter(sample_id == current_id) %>%
    st_buffer(buffer) %>%
    st_intersects(sf_obj) %>%
    unlist() %>%
    sf_obj[.,] %>%
    mutate(., nearest_feat = st_nearest_feature(.))
  current_feat <- neighbours %>%
    filter(sample_id == current_id)
  as.numeric(st_distance(current_feat, neighbours[current_feat$nearest_feat,]))
}
cl <- makeCluster(48)
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(dplyr))

# Calculate distance for n_samples
nn_dists <- map(n_samples, function(x){
  cat("n samples:", x, "\n")
  sub_sample <- sample_n(pixel_training_data, round(x * 0.2)) %>%
    select(sample_id, geometry)
  dists <- pblapply(sub_sample$sample_id,
                    dist_nearest_neighbour,
                    buffer = 10000, 
                    sf_obj = sub_sample, 
                    cl = cl) %>% 
    unlist()
  return(data.frame(n_samples = x,
                    mean_dist_nn = mean(dists, na.rm = T),
                    median_dist_nn = median(dists, na.rm = T),
                    min_dist_nn = min(dists, na.rm = T),
                    max_dist_nn = max(dists, na.rm = T),
                    n_na = sum(is.na(dists)))
  )
}) %>% 
  bind_rows()
nn_dists %>%
  gt() %>%
  gtsave("docs/nearest_neighbour_sensitivit.html")

# Calculate density for actual training dataset
rm(pixel_training_data)
load("data/training_data/pixel_training.Rda")
pixel_training_data <- pixel_training_data %>% 
  select(-contains("mask")) %>%
  select(-contains("250")) %>%
  select(-contains("mean_110")) %>%
  select(-contains("proportion")) %>%
  mutate(forest_value = factor(forest_value)) %>%
  na.omit()
st_crs(pixel_training_data)
pixel_training_data_high <- pixel_training_data %>% filter(forest_value == "high")
pixel_training_data_low <- pixel_training_data %>% filter(forest_value == "low")

dist_true_training_high <- pblapply(pixel_training_data_high$sample_id,
                               dist_nearest_neighbour,
                               buffer = 10000, 
                               sf_obj = pixel_training_data_high, 
                               cl = cl) %>% 
  unlist()
dist_true_training_low  <- pblapply(pixel_training_data_low$sample_id,
                                    dist_nearest_neighbour,
                                    buffer = 10000, 
                                    sf_obj = pixel_training_data_low, 
                                    cl = cl) %>% 
  unlist()
stopCluster(cl)
summary_stats_high <- data.frame(data_set = "pixel_training_data_high",
                                 mean_dist_nn = mean(dist_true_training_high, na.rm = T),
                                 median_dist_nn = median(dist_true_training_high, na.rm = T),
                                 min_dist_nn = min(dist_true_training_high, na.rm = T),
                                 max_dist_nn = max(dist_true_training_high, na.rm = T),
                                 n_na = sum(is.na(dist_true_training_high)))
                                 
summary_stats_low <- data.frame(data_set = "pixel_training_data_low",
                                 mean_dist_nn = mean(dist_true_training_low , na.rm = T),
                                 median_dist_nn = median(dist_true_training_low , na.rm = T),
                                 min_dist_nn = min(dist_true_training_low , na.rm = T),
                                 max_dist_nn = max(dist_true_training_low , na.rm = T),
                                 n_na = sum(is.na(dist_true_training_low )))

# > summary_stats_high
# data_set mean_dist_nn median_dist_nn min_dist_nn max_dist_nn n_na
# 1 pixel_training_data_high     94.59265       53.85165           0    8171.377    6
# > summary_stats_low
# data_set mean_dist_nn median_dist_nn min_dist_nn max_dist_nn n_na
# 1 pixel_training_data_low     77.43345       44.72136           0    9695.282    2