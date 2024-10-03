# Ranger model sensitivity to the sampling size

# Dependencies
library(caret)
library(tidyverse)
library(sf)
library(parallel)
library(pbapply)
library(ranger)
library(cowplot)
library(gt)

# Set target number of samples
# (adding an extra buffer of 30% more to allow for removal of NAs due to water)
n_samples_list <- c(paste0("_", c(10000, 45000, 60000, 75000) * 1.3), "")

# Prepare parallel envrionment
cl <- makeCluster(5)
clusterEvalQ(cl, library(ranger))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(caret))
clusterEvalQ(cl, library(ranger))

rf_fits <- pblapply(n_samples_list, cl = cl, function(n_samples){
  cat(n_samples, "\n")
  
  # Load data
  load(paste0("data/training_data/pixel_training_biowide", n_samples, ".Rda"))
  load(paste0("data/training_data/pixel_valid_biowide", n_samples, ".Rda"))
  
  # Set pseudo random generator seed
  set.seed(24231)

  # Prep data frames 
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
  
  # Set model hyperparemeter based on best preformance at 25 k
  hyper_params <- data.frame(mtry = 5,
                             min.node.size = 1,
                             splitrule	= "gini")
  
  # Fit ranodm forest and time it
  fit_time <- system.time(
    rf_fit <- train(forest_value ~ .,
                    data = train_data,
                    method = "ranger",
                    trControl = trainControl(method = "none",
                                             classProbs = TRUE, 
                                             summaryFunction = twoClassSummary),
                    tuneGrid = hyper_params,
                    importance = "permutation",
                    metric = "ROC",
                    num.threads = 1)
  )
  
  # Validate on test set
  test_preds <- predict(rf_fit, newdata = test_data)
  cf_fit <- confusionMatrix(data = test_preds, test_data$forest_value)
  
  # Return list with model, confusion matrix and summary data frame
  return(list(rf_fit, 
              cf_fit, 
              data.frame(
                n_samples = round((nrow(train_data) + nrow(test_data)) / 2),
                accuracy = cf_fit$overall[1],
                speed = fit_time[1]
              )))
})

stopCluster(cl)
rm(cl)

# Summarise results in overview table
(results_df <- bind_rows(map(rf_fits, function(x) x[[3]])) %>%
  arrange(n_samples) %>%
  remove_rownames())

# Export as table
results_df %>%
  mutate(accuracy = round(accuracy, 2),
         speed = round(speed)) %>%
  select(`Sapmples per class (n)` = n_samples, 
         Accuracy = accuracy,
         `Speed (s)` = speed) %>%
  gt() %>%
  tab_header(title = "random forest - Biowide") %>%
  gtsave("docs/sampling_size_sensitivity/sampling_size_sensitivity.html")

# Generate graph
(ranger_biowide_plot <- plot_grid(
  ggplot(results_df, aes(x = n_samples, y = accuracy)) +
    geom_line() +
    geom_point() +
    annotate("point", x = results_df$n_samples[2], y = results_df$accuracy[2],
             colour = "red") +
    annotate("text", x = 0, y = Inf, vjust = 1.5, hjust = 0,
             label = "random forest - biowide", size = 14 / .pt) +
    scale_x_continuous(limits = c(0, 100000)) +
    scale_y_continuous(limits = c(0.75, 0.9)) +
    labs(x = "Samples per class (n)", y = "Accuracy") +
    theme_cowplot(),
  ggplot(results_df, aes(x = n_samples, y = speed)) +
    geom_line() +
    geom_point() +
    annotate("point", x = results_df$n_samples[2], y = results_df$speed[2],
             colour = "red") +
    annotate("text", x = 0, y = Inf, vjust = 1.5, hjust = 0,
             label = "random forest - biowide", size = 14 / .pt) +
    scale_x_continuous(limits = c(0, 100000)) +
    scale_y_continuous(limits = c(0, 600)) +
    labs(x = "Samples per class (n)", y = "Training time (s)") +
    theme_cowplot()))
ranger_biowide_plot %>%
  save_plot("docs/sampling_size_sensitivity/sampling_sens_ranger_biowide.png",
            .,
            nrow = 1,
            ncol = 2,
            bg = "white")


## Calculate distances between samples

# helper function to calculate distances to nearest neighbour based on buffering
dist_nearest_neighbour <- function(current_id, buffer, sf_obj){
  # Retrieve current feature
  current_feat <- sf_obj %>%
    filter(sample_id == current_id)
  # Get all neighbours within buffer
  neighbours <-  sf_obj %>%
    filter(sample_id == current_id) %>%
    st_buffer(buffer) %>%
    st_intersects(sf_obj) %>%
    unlist() %>%
    sf_obj[.,] %>%
    filter(sample_id != current_id)
  # Identify nearest neighbour of current feature
  current_feat <- current_feat %>%
    mutate(nearest_feat = st_nearest_feature(., neighbours))
  # Calculate distance to nearest neighbour and return
  return(as.numeric(st_distance(current_feat, 
                                neighbours[current_feat$nearest_feat,])))
}

# Prepare parallel environment
cl <- makeCluster(48)
clusterEvalQ(cl, library(sf))
clusterEvalQ(cl, library(dplyr))

# Calculate distance for each set of n_samples
nn_dists <- map(n_samples_list, function(n_samples){
  cat("n samples:", n_samples, "\n")
  
  # Load data
  load(paste0("data/training_data/pixel_training", n_samples, ".Rda"))
  
  # Prep and clean as other data 
  pixel_training_data <- pixel_training_data  %>%
    select(-contains("mask")) %>%
    select(-contains("250")) %>%
    select(-contains("mean_110")) %>%
    select(-contains("proportion")) %>%
    select(-treetype_bjer_con) %>% 
    mutate(forest_value = factor(forest_value)) %>%
    na.omit()
  
  # Split into high and low 
  pixel_high <- filter(pixel_training_data, forest_value == "high")
  pixel_low <- filter(pixel_training_data, forest_value == "low")
  
  # Set pre-distance buffering
  buffer_dist <- 5000
  
  # Calculate distances to nearest neighbour within forest classes
  dists_high <- pblapply(pixel_high$sample_id,
                    dist_nearest_neighbour,
                    buffer = buffer_dist, 
                    sf_obj = pixel_high, 
                    cl = cl) %>% 
    unlist()
  dists_low <- pblapply(pixel_low$sample_id,
                         dist_nearest_neighbour,
                         buffer = buffer_dist, 
                         sf_obj = pixel_low, 
                         cl = cl) %>% 
    unlist()
  
  # Return summary statistics of ditsances
  return(data.frame(n_target = gsub("_", "", n_samples),
                    n_high = nrow(pixel_high),
                    n_low = nrow(pixel_low),
                    buffer_dist = buffer_dist,
                    mean_dist_high = mean(dists_high, na.rm = T),
                    median_dist_high = median(dists_high, na.rm = T),
                    min_dist_high = min(dists_high, na.rm = T),
                    max_dist_high = max(dists_high, na.rm = T),
                    n_na_high = sum(is.na(dists_high)),
                    mean_dist_low = mean(dists_low, na.rm = T),
                    median_dist_low = median(dists_low, na.rm = T),
                    min_dist_low = min(dists_low, na.rm = T),
                    max_dist_low = max(dists_low, na.rm = T),
                    n_na_low = sum(is.na(dists_low)))
  )
}) %>% 
  bind_rows() %>%
  arrange(n_high)

# Close parallel environment
stopCluster(cl)

# Export table to html
nn_dists %>%
  mutate(n_target = c(10000 * 1.3, 30000, 45000 * 1.3, 60000 * 1.3, 75000 * 1.3),
         prop_high = round(n_high / (n_low + n_high), 2)) %>%
  mutate(across(5:13, ~ round(.x))) %>%
  gt() %>%
  gtsave("docs/sampling_size_sensitivity/nearest_neighbour_distance.html")
