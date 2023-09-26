# Quick script to assess forest masks
# Jakob J. Assmann 24 August 2022 jakob.assmann@uzh.ch

library(terra)

# Bjerreskov forest mask
bjer_forest <- rast("data/predictor_data/treetype/forest_mask_bjer_above_half_ha.tif")

# Base map forest mask
basemap_forest <- rast("data/basemap_forests/forest_mask.tif")

# canopy height layer above 3m
canopy_height <- rast("F:/JakobAssmann/EcoDes-DK15_v1.1.0/canopy_height/canopy_height.vrt")
above_3m <- canopy_height >= 3
above_3m <- classify(above_3m, matrix(c(0, NA, 1, 1), ncol = 2, byrow = T))


# Sum up the area
area_above_3m <- length(cells(above_3m, 1)[[1]]) * 100
area_bjer_forest <- length(cells(bjer_forest, 1)[[1]]) * 100 
area_basemap_forest <- length(cells(basemap_forest, 1)[[1]]) * 100

# Add 3 m mask to bjer_forest mask
bjer_forest_3m <- mask(bjer_forest, above_3m)
writeRaster(bjer_forest_3m, "data/bjer_forest_3m.tif")
area_bjer_forest_3m <- length(cells(bjer_forest_3m, 1)[[1]]) * 100

# Apply minimum area mapping (500 m2)
forest_mask_polys <- as.polygons(bjer_forest_3m)
forest_mask_polys_sf <- st_as_sf(forest_mask_polys)
forest_mask_polys_sf <- st_cast(forest_mask_polys_sf, "POLYGON")
forest_mask_polys_sf_area <- st_area(forest_mask_polys_sf)
forest_mask_polys_sf_above_half_ha <- forest_mask_polys_sf[forest_mask_polys_sf_area >= units::set_units(500, "m^2"),]
forest_mask_bjer_3m_above_half_ha <- mask(bjer_forest_3m, vect(forest_mask_polys_sf_above_half_ha))

area_forest_mask_bjer_3m_above_half_ha <- length(cells(forest_mask_bjer_3m_above_half_ha, 1)[[1]]) * 100
-