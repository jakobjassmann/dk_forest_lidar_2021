# Comparsion between the two forest type variables
# Jakob J. Assmann 16 February 2022

# Housekeeping
library(tidyverse)
library(ggplot2)
library(sf)

# Load data
load("data/training_data/pixel_training.Rda")

## Compare broadleaf forests

# Confusion matrix
pixel_training_data %>%
  st_drop_geometry() %>%
  select(forest_type_dec, treetype_bjer_dec) %>%
  table()
#             treetype_bjer_dec
# forest_type_dec      0      1
#               0  27635   5575
#               1  11873 119273


# Accuracy
sum(pixel_training_data$forest_type_dec == pixel_training_data$treetype_bjer_dec, na.rm = T) /
   nrow(pixel_training_data)
# [1] 0.73454

# Sensitivity (treetype_bjer = true values)
pixel_training_data %>% 
  st_drop_geometry() %>%
  select(forest_type_dec, treetype_bjer_dec) %>%
  filter(treetype_bjer_dec == 1) %>%
  na.omit() %>%
  summarize(sens = sum(forest_type_dec == treetype_bjer_dec) / n())
#   sens
# 1 0.9553457

# Specificity (treetype_bjer = true values)
pixel_training_data %>% 
  st_drop_geometry() %>%
  select(forest_type_dec, treetype_bjer_dec) %>%
  filter(treetype_bjer_dec == 0) %>%
  na.omit() %>%
  summarize(spec = sum(forest_type_dec == treetype_bjer_dec) / n())
# sens
# 1 0.6994786

## Compare confierous forests

# Confusion matrix
pixel_training_data %>%
  st_drop_geometry() %>%
  select(forest_type_con, treetype_bjer_con) %>%
  table()
#           treetype_bjer_con
# forest_type_con      0      1
#               0 119435  12223
#               1   5413  27285

# Accuracy
sum(pixel_training_data$forest_type_con == pixel_training_data$treetype_bjer_con, na.rm = T) /
  nrow(pixel_training_data)
# [1] 0.73454

# Sensitivity (treetype_bjer = true values)
pixel_training_data %>% 
  st_drop_geometry() %>%
  select(forest_type_con, treetype_bjer_con) %>%
  filter(treetype_bjer_con == 1) %>%
  na.omit() %>%
  summarize(sens = sum(forest_type_con == treetype_bjer_con) / n())
#   sens
# 1 0.6906196

# Specificity (treetype_bjer = true values)
pixel_training_data %>% 
  st_drop_geometry() %>%
  select(forest_type_con, treetype_bjer_con) %>%
  filter(treetype_bjer_con == 0) %>%
  na.omit() %>%
  summarize(spec = sum(forest_type_con == treetype_bjer_con) / n())
#   spec
# 1 0.9566433