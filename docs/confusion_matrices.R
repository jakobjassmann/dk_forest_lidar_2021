# Generate confusion matrix tables
# Jakob J. Assmann jakob.assmann@uzh.ch 4 June 2024

# Dependencies
library(tidyverse)
library(gt)
library(caret)
library(ranger)
library(gbm)

# Load models
load("data/models/final_ranger_model_pixel_biowide.Rda")
rf_biowide <- rf_fit
rm(rf_fit)
load("data/models/final_ranger_model_pixel_derek.Rda")
rf_derek <- rf_fit
rm(rf_fit)
load("data/models/final_gbm_model_pixel_biowide.Rda")
gbm_biowide <- gbm_fit
rm(gbm_fit)
load("data/models/final_gbm_model_pixel_derek.Rda")
gbm_derek <- gbm_fit
rm(gbm_fit)

# Load validation data
load("data/training_data/pixel_valid_biowide.Rda")
load("data/training_data/pixel_valid_derek.Rda")

# Get confusion matrices
cf_rf_biowide <- confusionMatrix(predict(rf_biowide, 
                                         newdata = pixel_valid_biowide),
                pixel_valid_biowide$forest_value)
cf_rf_derek <- confusionMatrix(predict(rf_derek, 
                                         newdata = pixel_valid_derek),
                                 pixel_valid_derek$forest_value)
cf_gbm_biowide <- confusionMatrix(predict(gbm_biowide, 
                                         newdata = pixel_valid_biowide),
                                 pixel_valid_biowide$forest_value)
cf_gbm_derek <- confusionMatrix(predict(gbm_derek, 
                                         newdata = pixel_valid_derek),
                                 pixel_valid_derek$forest_value)

# Generate confusion matrix table for the ranger models
list(cbind(cf_rf_biowide$table[,1], cf_rf_biowide$table[,2]) %>%
       as.data.frame() %>%
       remove_rownames() %>%
       mutate(prediction = c("high", "low")) %>%
       relocate(prediction) %>%
       set_names(c("Prediction", "high", "low")) %>%
       gt(rowname_col = "Prediction") %>%
       tab_stubhead("Prediction") %>%
       tab_spanner("Reference", columns = c("high", "low")) %>%
       as_raw_html(),
     cbind(cf_rf_derek$table[,1], cf_rf_derek$table[,2]) %>%
       as.data.frame() %>%
       remove_rownames() %>%
       mutate(prediction = c("high", "low")) %>%
       relocate(prediction) %>%
       set_names(c("Prediction", "high", "low")) %>%
       gt(rowname_col = "Prediction") %>%
       tab_stubhead("Prediction") %>%
       tab_spanner("Reference", columns = c("high", "low")) %>%
       as_raw_html()) %>%
  data.frame() %>%
  setNames(c("Random Forest - Biowide",
             "Random Forest - SustainScapes")) %>%
  gt() %>%
  fmt_markdown(columns = everything()) %>%
  gtsave("docs/confusion_matrices/ranger_models.html")

# Generate confusion matrix table for the gbm models
list(cbind(cf_gbm_biowide$table[,1], cf_gbm_biowide$table[,2]) %>%
       as.data.frame() %>%
       remove_rownames() %>%
       mutate(prediction = c("high", "low")) %>%
       relocate(prediction) %>%
       set_names(c("Prediction", "high", "low")) %>%
       gt(rowname_col = "Prediction") %>%
       tab_stubhead("Prediction") %>%
       tab_spanner("Reference", columns = c("high", "low")) %>%
       as_raw_html(),
     cbind(cf_gbm_derek$table[,1], cf_gbm_derek$table[,2]) %>%
       as.data.frame() %>%
       remove_rownames() %>%
       mutate(prediction = c("high", "low")) %>%
       relocate(prediction) %>%
       set_names(c("Prediction", "high", "low")) %>%
       gt(rowname_col = "Prediction") %>%
       tab_stubhead("Prediction") %>%
       tab_spanner("Reference", columns = c("high", "low")) %>%
       as_raw_html()) %>%
  data.frame() %>%
  setNames(c("gbm - Biowide",
             "gbm - SustainScapes")) %>%
  gt() %>%
  fmt_markdown(columns = everything()) %>%
  gtsave("docs/confusion_matrices/gbm_models.html")

