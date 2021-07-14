# Preparation script to extract EcoDes-DK15 data for the 2016 Basemap polygons.
# Jakob J. Assmann 13 July 2021

# Dependencies
library(raster)
library(exactextractr)
library(parallel)
library(sf)
library(tidyverse)

# Load basemap forest polygons
forest_polys <- read_sf("data/projections/basemap2016_forest.shp")

# Load list of Ecodes-DK variabels
ecodes_vrt <- read.table("D:/Jakob/dk_nationwide_lidar/data/outputs/list_of_vrts.txt",
                         stringsAsFactors = F)[,1]
# Remove unneded variables
ecodes_vrt <- ecodes_vrt[!grepl("point_count",ecodes_vrt)]
ecodes_vrt <- ecodes_vrt[!grepl("point_source",ecodes_vrt)]
ecodes_vrt <- ecodes_vrt[!grepl("building",ecodes_vrt)]

# add full file name
ecodes_vrt <- paste0("D:/Jakob/dk_nationwide_lidar/data/outputs/", ecodes_vrt)

# Transform polygons to raster data CRS
forest_polys <- st_transform(forest_polys,
                               st_crs(raster("D:/Jakob/dk_nationwide_lidar/data/outputs/dtm_10m/dtm_10m_6049_684.tif")))

# Extract mean statistics for polygons using exactextract
cl <- makeCluster(42)
clusterEvalQ(cl, library(exactextractr))
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(sf))
system.time(poly_projection_data <- parLapply(cl, 
                                            ecodes_vrt, 
                                            function(vrt_file, sample_polys){
                                              ecodes_raster <- raster(vrt_file)
                                              sample_polys$mean <- exact_extract(ecodes_raster, sample_polys, fun = "mean")
                                              sample_polys$sd <- exact_extract(ecodes_raster, sample_polys, fun = "stdev")
                                              names(sample_polys)[4] <- paste0(gsub(".*/(.*).vrt", "\\1", vrt_file), "_mean")
                                              names(sample_polys)[5] <- paste0(gsub(".*/(.*).vrt", "\\1", vrt_file), "_sd")
                                              return(sample_polys)
                                            },
                                            forest_polys))
save(poly_projection_data, file = "data/projections/poly_projections_data.Rda")
stopCluster(cl)

# Combine dataframes
forest_polys <-  mutate(forest_polys, forest_id = object_id) %>% select(-DN, -object_id)
poly_projection_data_tidy <- poly_projection_data %>% map(st_drop_geometry) %>%
  map(function(x) mutate(x, forest_id = object_id) %>% select(-DN, -object_id)) %>% 
  reduce(full_join, by = "forest_id") %>% full_join(forest_polys, ., by = "forest_id") 

# Safe tidy data
save(poly_projection_data_tidy, file = "data/projections/poly_projections_data_tidy.Rda")
