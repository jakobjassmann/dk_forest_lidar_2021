# Reporject predictor data where required
# Jakob J. Assmann j.assmann@bio.au.dk 16 February 2022

# Housekeeping
library(terra)
terraOptions(progress = 1)

# Load rasters
clay <- rast("data/predictor_data/soil_layers/Clay.tif")
sand <- rast("data/predictor_data/soil_layers/Sand.tif")
soil_carbon <- rast("data/predictor_data/soil_layers/Soc.tif")
ns_groundwater_summer <- rast("data/predictor_data/terraennaert_grundvand_10m/summer_predict.tif")

# Load target raster
dtm10m <- rast("F:/JakobAssmann/EcoDes-DK_v1.1.0/dtm_10m/dtm_10m.vrt")

# Reproject
project(clay, dtm10m, filename = "data/predictor_data/soil_layers/Clay_utm32_10m.tif",
        method = "near")
project(sand, dtm10m, filename = "data/predictor_data/soil_layers/Sand_utm32_10m.tif",
        method = "near")
project(soil_carbon, dtm10m, filename = "data/predictor_data/soil_layers/Soc_utm32_10m.tif",
        method = "near")

# Set CRS for ground water raster and reproject also 
crs(ns_groundwater_summer) <- "EPSG:25832"
project(ns_groundwater_summer, dtm10m, 
        filename ="data/predictor_data/terraennaert_grundvand_10m/ns_groundwater_summer_utm32_10m.tif",
        method = "bilinear")
