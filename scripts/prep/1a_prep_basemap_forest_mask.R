# Quick script to prepare basemap forest mask
# Jakob J. Assmann j.assmann@bio.au.d 22 February 2022

# Dependencies
library(terra)

# Load forest mask
sub_tree_cover <- rast("data/basemap_forests/sub_tree_cover_2016.tif")
forest_mask <- sub_tree_cover == 2

# Project to EcoDes-DK15 grid
dtm10m <- rast("F:/JakobAssmann/EcoDes-DK15_v1.1.0/dtm_10m/dtm_10m.vrt")
forest_mask <- terra::project(forest_mask, dtm10m, method = "near")
forest_mask <- crop(forest_mask, dtm10m)

# Set 0 to NA
forest_mask <- classify(forest_mask, 
                        matrix(c(0, NA, 
                                 1, 1),
                               byrow = T, ncol = 2))
writeRaster(forest_mask,
            filename = "data/basemap_forests/forest_mask.tif")
