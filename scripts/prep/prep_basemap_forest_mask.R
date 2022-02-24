# Quick script to prepare basemap forest mask
# Jakob J. Assmann j.assmann@bio.au.d 22 February 2022

# Dependencies
library(terra)

# Load forest mask
sub_tree_cover <- rast("data/basemap_forests/sub_tree_cover_2016.tif")
forest_mask <- sub_tree_cover == 2

dtm10m <- rast("F:/JakobAssmann/EcoDes-DK15_v1.1.0/dtm_10m/dtm_10m.vrt")
forest_mask <- project(forest_mask, dtm10m, method = "near")
writeRaster(forest_mask,
            filename = "data/basemap_forests/forest_mask.tif")
