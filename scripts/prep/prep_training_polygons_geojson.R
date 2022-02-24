# Quick script to convert training polygons to json files for leaflet
# Jakob J. Assmann j.assmann@bio.au.dk 24 February 2022

# Load dependencies
library(tidyverse)
library(sf)

# Load polygons
load("data/training_data/polygon_geometries/high_quality_polys.Rda")
load("data/training_data/polygon_geometries/low_quality_polys.Rda")

# Set properties
high_quality$forest_quality <- "high"
low_quality$forest_quality <- "low"
combined <- bind_rows(high_quality, low_quality) %>%
  st_transform(crs = 4326)
# Export
write_sf(combined, dsn = "data/projections/training_polygons.geojson")
write_sf(combined[1:10,],
         dsn = "data/projections/training_polygons_small.geojson")
