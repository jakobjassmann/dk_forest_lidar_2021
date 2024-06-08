# Stratified split of training data into training and validation for the
# sample size sensitivity analysis

# Jakob J. Assmann jakob.assmann@uzh.ch 1 June 2024

## 1) Preparations ----

# Dependencies
library(tidyverse)
library(sf)

# Set target number of samples
# (adding an extra buffer of 30% more to allow for removal of NAs due to water)
n_samples <- c(10000, 45000, 60000, 75000) * 1.3

# map over sample sizes
map(n_samples, function(n_samples){
  # Load data
  load(paste0("data/training_data/pixel_training_", n_samples, ".Rda"))
  
  # Filter out masked data and na values   
  pixel_training_data <- pixel_training_data %>% 
    st_drop_geometry() %>%
    select(-contains("mask")) %>%
    select(-contains("250")) %>%
    select(-contains("mean_110")) %>%
    select(-contains("proportion")) %>% 
    select(-treetype_bjer_con) %>%
    mutate(forest_value = factor(forest_value)) %>%
    na.omit()
  
  # Check VIFs (takes a while to correlate)
  usdm::vif(pixel_training_data[,-(1:4)])
  
  # Set seed for pseudo random numbers
  set.seed(70222)
  
  ## 2) Biowide stratification ----
  
  # Group data
  pixel_training_data_biowide <- group_by(pixel_training_data,
                                          biowide_region) 
  
  # Random strafiied Split 80% training, 20% validation
  pixel_training_biowide <-  sample_frac(pixel_training_data_biowide,
                                         size = 0.8)
  pixel_valid_biowide <- pixel_training_biowide %>%
    ungroup() %>%
    pull(sample_id) %>%
    as.character() %>%
    {as.character(pixel_training_data_biowide$sample_id) %in% .} %>%
    `!` %>%
    pixel_training_data_biowide[.,]
  
  # Save files
  save(pixel_training_biowide, 
       file = paste0("data/training_data/pixel_training_biowide_", n_samples, ".Rda"))
  save(pixel_valid_biowide,
       file = paste0("data/training_data/pixel_valid_biowide_", n_samples, ".Rda"))
  
  ## 2) Derek's stratification ----
  
  # Group data
  pixel_training_data_derek <- group_by(pixel_training_data, 
                                        dereks_stratification) 
  
  # Random strafiied Split 80% training, 20% validation
  pixel_training_derek <-  sample_frac(pixel_training_data_derek,
                                       size = 0.8)
  pixel_valid_derek <- pixel_training_derek %>%
    ungroup() %>%
    pull(sample_id) %>%
    as.character() %>%
    {as.character(pixel_training_data_derek$sample_id) %in% .} %>%
    `!` %>%
    pixel_training_data_derek[.,]
  
  # Save files
  save(pixel_training_derek, 
       file = paste0("data/training_data/pixel_training_derek_", n_samples, ".Rda"))
  save(pixel_valid_derek,
       file = paste0("data/training_data/pixel_valid_derek_", n_samples, ".Rda"))
})
