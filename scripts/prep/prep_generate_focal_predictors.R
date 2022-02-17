# DK Forest LiDAR - Generate focal predictor raster 
# Jakob J. Assmann j.assmann@bio.au.dk 16 February 2022

# Dependencies
library(terra)
library(pbapply)
library(parallel)

# Load rasters
dtm_10m  <- rast("F:/JakobAssmann/EcoDes-DK_v1.1.0/dtm_10m/dtm_10m.vrt")
ns_groundwater_summer <- rast("data/predictor_data/terraennaert_grundvand_10m/summer_predict.tif")
canopy_height <- rast("F:/JakobAssmann/EcoDes-DK_v1.1.0/canopy_height/canopy_height.vrt")
vegetation_density <-  rast("F:/JakobAssmann/EcoDes-DK_v1.1.0/vegetation_density/vegetation_density.vrt")

# Prepare parameters for focal call
focal_params <- data.frame(
  raster_file = c(rep("F:/JakobAssmann/EcoDes-DK_v1.1.0/dtm_10m/dtm_10m.vrt", 4),
                  rep("data/predictor_data/terraennaert_grundvand_10m/ns_groundwater_summer_utm32_10m.tif", 4),
                  rep("F:/JakobAssmann/EcoDes-DK_v1.1.0/canopy_height/canopy_height.vrt", 4),
                  rep("F:/JakobAssmann/EcoDes-DK_v1.1.0/vegetation_density/vegetation_density.vrt", 4)),
  raster_name = c(rep("dtm_10m", 4),
                  rep("ns_groundwater_summer", 4),
                  rep("canopy_height", 4),
                  rep("vegetation_density", 4)),
  window = rep(c(11,25), 8),
  fun = rep(c("mean", "mean", "sd", "sd"), 4)
)
focal_params$out_file <- paste0("data/predictor_data/focal_variables/",
                                focal_params$raster_name,
                                "_", focal_params$fun, 
                                "_", focal_params$window * 10, "m.tif")


# Calculate focal parameters in parallel
cl <- makeCluster(16)
clusterEvalQ(cl, library(terra))

pblapply(
  split(focal_params, paste0(focal_params$raster_name,
    "_", focal_params$fun, 
    "_", focal_params$window * 10, "m")),
  function(params){
    focal(rast(params$raster_file),
          w = params$window,
          fun = params$fun,
          na.policy = "omit",
          filename = params$out_file,
          overwrite = T)
    return(0)
    },
  cl = cl)
