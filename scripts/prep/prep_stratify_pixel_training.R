# Stratified split of training data into training and validation
# Jakob J. Assmann j.assmann@bio.au.dk

## 1) Preparations ----

# Dependencies
library(tidyverse)
library(sf)

# Load data
load("data/training_data/pixel_training.Rda")

# Filter out masked data and na values   
pixel_training_data <- pixel_training_data %>% 
  st_drop_geometry() %>%
  filter(!is.na(inland_water_mask)) %>%
  filter(!is.na(sea_mask)) %>%
  select(-contains("mask")) %>%
  select(-contains("250")) %>%
  select(-contains("mean_110")) %>%
  select(-contains("proportion")) %>%
  select(-heat_load_index,
         -aspect,
         -openness_mean,
         -normalized_z_mean,
         -canopy_openness,
         -treetype_bjer_con,
         -twi) %>% 
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
     file = "data/training_data/pixel_training_biowide.Rda")
save(pixel_valid_biowide,
     file = "data/training_data/pixel_valid_biowide.Rda")

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
     file = "data/training_data/pixel_training_derek.Rda")
save(pixel_valid_derek,
     file = "data/training_data/pixel_valid_derek.Rda")
