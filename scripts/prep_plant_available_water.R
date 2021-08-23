# Prep script to get the plant avaialble water variable ready for the predictions
# Jakob J. Assmann j.assmann@bio.au.dk 13 August 2021

# Dependencies
library(terra)
library(tidyverse)

# Load files
paw_stack <- list.files("data/plant_available_water/", pattern = ".tif$", full.names = T) %>%
  map(rast) %>% rast()

# Calculate total paw for 160 cm
paw_160cm <- 3 * paw_stack[[1]] + 3 * paw_stack[[2]] + 4* paw_stack[[3]] + 6* paw_stack[[4]]
writeRaster(paw_160cm, "data/plant_available_water/paw_160cm.tif")
