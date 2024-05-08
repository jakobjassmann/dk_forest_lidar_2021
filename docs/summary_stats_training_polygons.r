# Sumarry stats of training data for supplementary material
# Jakob J. Assmann jakob.assmann@uzh.ch 30 November 2023
library(tidyverse)
library(sf)

load("data/training_data/polygon_geometries/high_quality_polys.Rda")
load("data/training_data/polygon_geometries/low_quality_polys.Rda")

high_quality %>%
    mutate(area = st_area(geometry)) %>%
    st_drop_geometry() %>%
    group_by(polygon_source) %>%
    summarize(mean(area) *10^(-6), sd(area)*10^(-6))

low_quality %>%
    mutate(area = st_area(geometry)) %>%
    st_drop_geometry() %>%
    group_by(polygon_source) %>%
    summarize(mean(area) *10^(-6), sd(area)*10^(-6))
0.010