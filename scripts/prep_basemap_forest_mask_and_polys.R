# Wee script to prepare forest basemap for bornholm.
library(raster)
library(rgdal)
library(sf)
library(fasterize)
library(parallel)

# Tweak raster environment
rasterOptions(progress = "text")
rasterOptions(maxmemory = 5e+11)
rasterOptions()

# Load basemap 
base_map <- raster("data/projections/Basemap02_public/export/LU_02_20161.tif")


# turn the two forest layers (110000 and 110110) into a binary mask
base_map <- reclassify(base_map, c(1, 109999, NA, 
                                   109999, 110001, 1,
                                   110001, 110109, NA,
                                   110109, 110111, 1,  
                                   110111, 999000, NA
                                   ))

NAvalue(base_map) <- 0

# export raster
writeRaster(base_map,
            "data/projections/basemap2016_forest_mask.tif",
            overwrite = T)
# 
# # Load shapefiles for p15 and p25 forests
# p15_shp <- shapefile("data/p15_forests/shapefiles/skov_kortlaegning_2016_2018.shp")
# p25_shp <- shapefile("data/p25_forests/shapefiles/p25_offentligareal.shp")
# p15_shp <- spTransform(p15_shp, crs(base_map))
# p25_shp <- spTransform(p25_shp, crs(base_map))
# p15_shp$object_id <- paste0("p15_", 1:length(p15_shp))
# p25_shp$object_id <- paste0("p25_", 1:length(p25_shp))
# p15_shp <- p15_shp[,"object_id"]
# p25_shp <- p25_shp[,"object_id"]
# 
# # Combine and devide into chuncks of 54
# p_combined <- rbind(p15_shp, p25_shp)
# p_combined <- st_as_sf(p_combined)
# base_mask <- fasterize(p_combined, base_map)
# writeRaster(base_mask, "data/nationwide_forests/Basemap02_public/export/base_mask.tif")
# 
# # mask p15 and p25 forests form base map
# base_map <- mask(base_map, base_mask, inverse = T)
# writeRaster(base_map, "data/nationwide_forests/Basemap02_public/export/basemap_forest_raster.tif",
#             overwrite = T)

# Polygonize raster using gdal through OSGeo4W64 (this is way faster than doing it in R)
# CAREFUL! This will likely take a while - around 1h or so.
system2(command = "C:/OSGeo4W64/OSGeo4W.bat",
        args = c("gdal_polygonize",
                 paste0(getwd(), "/data/projections/basemap2016_forest_mask.tif"),
                 "-b 1",  "-f \"ESRI Shapefile\"",
                 paste0(getwd(), "/data/projections/basemap2016_forest.shp") 
                 ))

# Load as shapefile
dk_basemap_forests <- shapefile(paste0(getwd(), "/data/projections/basemap2016_forest.shp"))

# add object ID
dk_basemap_forests$object_id <- paste0("basemap_" ,
                                     1:length(dk_basemap_forests@polygons))

# # Prep cluster
# cl <- makeCluster(54)
# 
# ## Load packages into cluster
# clusterEvalQ(cl, library(raster))
# clusterEvalQ(cl, library(sp))
# clusterExport(cl, "dk_base_map_forests")
# # filter by minmum area (300 m2)
# polys_to_keep <- parLapply(cl, dk_base_map_forests$object_id, function(x){
#   poly <- dk_base_map_forests[which(dk_base_map_forests$object_id == x),]
#   if (area(poly) < 300) {
#     return(F)
#   } else return(T)
# })
# 
# polys_to_keep <- unlist(polys_to_keep)
# dk_base_map_forests <- dk_base_map_forests[polys_to_keep,] 
# 
# stopCluster(cl)
# 
# # re-assign object ID
# dk_base_map_forests$object_id <- paste0("DK_forest_" ,
#                                      1:length(dk_base_map_forests@polygons))

# export
shapefile(dk_basemap_forests, 
          paste0(getwd(), "/data/projections/basemap2016_forest.shp"),
          overwrite = T)
