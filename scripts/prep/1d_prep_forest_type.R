# Helper script to prepare foresty type rasters for predictions
# Jakob J Assmann j.assmann@bio.au.dk 12 August 2021

# # Dependencies
# library(raster)
# library(sf)

# # Load forest type raster for mainland DK
# forest_type <- raster("data/predictor_data/conif_vs_broadleaf/dk_forest_con_vs_dec.tif")
# 
# # Create factor dummy rasters
# writeRaster(forest_type == 1, "data/conif_vs_broadleaf/forest_type_cloud.tif", 
#             overwrite = T)
# writeRaster(forest_type == 2, "data/conif_vs_broadleaf/forest_type_con.tif", 
#             overwrite = T)
# writeRaster(forest_type == 3, "data/conif_vs_broadleaf/forest_type_dec.tif", 
#             overwrite = T)
# 
# # Prepare and add Bornholm forest type classification
# bornholm <- raster("data/predictor_data/conif_vs_broadleaf/bornholm_forest_class.tif")
# bornholm[bornholm[] == 2] <- 3
# bornholm[bornholm[] == 1] <- 2
# bornholm <- ratify(bornholm)
# levels(bornholm) <- levels(forest_type)
# bornholm_extent_utm32 <- projectExtent(bornholm, forest_type)
# bornholm <- projectRaster(bornholm, bornholm_extent_utm32, method = "ngb")
# plot(bornholm)
# writeRaster(bornholm, "data/predictor_data/conif_vs_broadleaf/bornholm_forest_con_vs_dec.tif", overwrite = T)
# writeRaster(bornholm == 1, "data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_cloud.tif", overwrite = T)
# writeRaster(bornholm == 2, "data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_con.tif", overwrite = T)
# writeRaster(bornholm == 3, "data/predictor_data/conif_vs_broadleaf/bornholm_forest_type_dec.tif", overwrite = T)

# Bjerreskov et al. 2021 forest types
library(terra)
terraOptions(progress = 1)

# Load raster and create binary masks
treetype_bjer <- rast("data/predictor_data/treetype/treetype_dk_raw.tif")
treetype_bjer_con <- treetype_bjer == 1
treetype_bjer_dec <- treetype_bjer == 2

# Load target raster for resampling
dtm10m <- rast("F:/JakobAssmann/EcoDes-DK15_v1.1.0/dtm_10m/dtm_10m.vrt")

# Resample
treetype_bjer_con <- resample(treetype_bjer_con, dtm10m, method = "near")
treetype_bjer_dec <- resample(treetype_bjer_dec, dtm10m, method = "near")

# Write out rasters
writeRaster(treetype_bjer_con, "data/predictor_data/treetype/treetype_bjer_con.tif")
writeRaster(treetype_bjer_dec, "data/predictor_data/treetype/treetype_bjer_dec.tif")

# Generate alternative forest mask
forest_mask <- resample(!is.na(treetype_bjer), dtm10m, method = "near")
forest_mask <- classify(forest_mask, 
                        matrix(c(0, NA, 
                                 1, 1),
                               byrow = T, ncol = 2))

# Apply minimum area mapping (500 m2)
forest_mask_polys <- as.polygons(forest_mask)
forest_mask_polys_sf <- st_as_sf(forest_mask_polys)
forest_mask_polys_sf <- st_cast(forest_mask_polys_sf, "POLYGON")
forest_mask_polys_sf_area <- st_area(forest_mask_polys_sf)
forest_mask_polys_sf_above_half_ha <- forest_mask_polys_sf[forest_mask_polys_sf_area >= units::set_units(500, "m^2"),]
forest_mask_above_half_ha <- mask(forest_mask, vect(forest_mask_polys_sf_above_half_ha))

# Write forest masks
writeRaster(forest_mask,
            filename = "data/predictor_data/treetype/forest_mask_bjer.tif")
writeRaster(forest_mask_above_half_ha,
            filename = "data/predictor_data/treetype/forest_mask_bjer_above_half_ha.tif")

# EOF