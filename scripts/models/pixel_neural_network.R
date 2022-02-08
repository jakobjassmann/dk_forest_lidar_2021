# First attempt to train a neural network to predict the forest class
# Jakob Assmann j.assmann@bio.au.dk 8 July 2021

# Dependencies
library(keras)
library(tidyverse)

# Load data
load("data/pixel_sample.Rda")

# Reshape data 
combined_sample <- combined_sample %>% map(function(x) dplyr::select(x, -forest_class)) %>%
  reduce(full_join, by = "sample_id") %>% full_join(combined_sample[[1]][,1:2], 
                                                    ., by = "sample_id")
# Filter out masked data and na values   
combined_sample <- combined_sample %>% 
  filter(!is.na(inland_water_mask)) %>%
  filter(!is.na(sea_mask)) %>%
  select(-contains("mask")) %>%
  mutate(forest_class = factor(forest_class)) %>%
  na.omit()

# Split data and convert to matrix
one_fifth <- round(nrow(combined_sample)*.2)
one_tenth <- round(nrow(combined_sample)*.1)
valid_index <- sample(1:nrow(combined_sample), one_tenth, replace = F)
test_index <- sample((1:nrow(combined_sample))[-valid_index], one_fifth, replace = F)

train_data <- combined_sample[c(-valid_index, -test_index),c(-1,-2)] %>% as.matrix()
valid_data <- combined_sample[valid_index,c(-1,-2)] %>% as.matrix()
test_data <- combined_sample[test_index,c(-1,-2)] %>% as.matrix()

train_labels <- combined_sample[c(-valid_index, -test_index),2] %>% as.numeric() -1
valid_labels <- combined_sample[valid_index,2] %>% as.numeric() -1
test_labels <- combined_sample[test_index,2] %>% as.numeric() -1

# Standardise predictors based on training data
means <- colMeans(train_data)
stds <- apply(train_data, 2, sd)

train_data <- scale(train_data, center = means, scale = stds)
valid_data <- scale(valid_data, center = means, scale = stds)
test_data <- scale(test_data, center = means, scale = stds)

# Define Network
network <- keras_model_sequential() %>%
  layer_dense(units = 512, activation = "relu", input_shape = 39) %>%
  layer_dense(units = 256, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")

# Check network
summary(network)

# Compile network
network %>% compile(
  optimizer = "rmsprop",
  loss = "binary_crossentropy",
  metrics = c("accuracy")
)

# Train network for 20 epochs to start with
history <- network %>% fit(
  train_data,
  train_labels,
  epochs = 20,
  batch_size = 256,
  validation_data = list(valid_data, valid_labels)
)

# Fit final network with 6 epochs
# Define Network
network <- keras_model_sequential() %>%
  layer_dense(units = 512, activation = "relu", input_shape = 39) %>%
  layer_dense(units = 256, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")

# Compile network
network %>% compile(
  optimizer = "rmsprop",
  loss = "binary_crossentropy",
  metrics = c("accuracy")
)

# Train network for 20 epochs to start with
history <- network %>% fit(
  train_data,
  train_labels,
  epochs = 6,
  batch_size = 256,
  validation_data = list(valid_data, valid_labels)
)

# Validate
predictions <- network %>% predict_classes(test_data)
mean(predictions == test_labels)                    
(conf_matrix <- table(test_labels, predictions))
(sens <- conf_matrix[1,1] / sum(conf_matrix[1,]))
(spec <- conf_matrix[2,2] / sum(conf_matrix[2,]))
(tss <- sens + spec - 1)
