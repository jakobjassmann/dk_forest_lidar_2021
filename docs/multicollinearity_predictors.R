library(tidyverse)
library(sf)
load("data/training_data/pixel_training_biowide.Rda")
training_vifs <- pixel_training_biowide %>%
    ungroup() %>%
    st_drop_geometry() %>%
    select(-(1:6)) %>%
    data.frame() %>%
    usdm::vif()
training_cors <- pixel_training_biowide %>%
    ungroup() %>%
    st_drop_geometry() %>%
    select(-(1:6)) %>%
    data.frame() %>%
    cor()
