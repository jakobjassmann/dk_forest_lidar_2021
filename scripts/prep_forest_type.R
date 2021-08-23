# Helper script to prepare foresty type rasters for predictions
# Jakob J Assmann j.assmann@bio.au.dk 12 August 2021

# Dependencies
library(raster)
library(sf)

# Load forest type raster for mainland DK
forest_type <- raster("data/conif_vs_broadleaf/dk_forest_con_vs_dec.tif")

# Create factor dummy rasters
writeRaster(forest_type == 1, "data/conif_vs_broadleaf/forest_type_cloud.tif", 
            overwrite = T)
writeRaster(forest_type == 2, "data/conif_vs_broadleaf/forest_type_con.tif", 
            overwrite = T)
writeRaster(forest_type == 3, "data/conif_vs_broadleaf/forest_type_dec.tif", 
            overwrite = T)

# Prepare and add Bornholm forest type classification
bornholm <- raster("data/conif_vs_broadleaf/bornholm_forest_class.tif")
bornholm[bornholm[] == 2] <- 3
bornholm[bornholm[] == 1] <- 2
bornholm <- ratify(bornholm)
levels(bornholm) <- levels(forest_type)
bornholm_extent_utm32 <- projectExtent(bornholm, forest_type)
bornholm <- projectRaster(bornholm, bornholm_extent_utm32)
plot(bornholm)
writeRaster(bornholm, "data/conif_vs_broadleaf/bornholm_forest_con_vs_dec.tif", overwrite = T)
writeRaster(bornholm == 1, "data/conif_vs_broadleaf/bornholm_forest_type_cloud.tif")
writeRaster(bornholm == 2, "data/conif_vs_broadleaf/bornholm_forest_type_con.tif")
writeRaster(bornholm == 3, "data/conif_vs_broadleaf/bornholm_forest_type_dec.tif")


