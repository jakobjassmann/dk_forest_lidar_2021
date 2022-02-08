# Script to prepare neighbourhood (focal) variables for set of most important variables
# Jakob Assmann j.assmann@bio.au.dk 24 August 2021

# Depenencies
library(terra)
library(parallel)

# Load variables to be computed
normalized_z_sd <- "D:/Jakob/dk_nationwide_lidar/data/outputs/normalized_z_sd/normalized_z_sd.vrt"
paw_160cm <- "data/plant_available_water/paw_160cm.tif"
canopy_height <- "D:/Jakob/dk_nationwide_lidar/data/outputs/canopy_height/canopy_height.vrt"

# Set up cluster
cl <- makeCluster(3)
clusterEvalQ(cl, library(terra))

# Run focal analysis
parLapply(cl,
          list(normalized_z_sd,
               paw_160cm,
               canopy_height),
          function(rast_file){
            rast_object <- rast(rast_file)
            stat_fun <- "sd"
            if(names(rast_object) == "normalized_z_sd ") stat_fun <- "mean"
            focal_rast <- focal(rast_object, w=3, fun=stat_fun, na.rm=FALSE)
            writeRaster(focal_rast, paste0("data/focal_variables/", names(rast_object), "_focal_3x3.tif"))
            return("Done.")
          })

# Stop cluster
stopCluster(cl)