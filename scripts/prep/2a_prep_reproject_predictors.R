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
dtm10m <- rast("F:/JakobAssmann/EcoDes-DK15_v1.1.0/dtm_10m/dtm_10m.vrt")

# Expand / gap fill soil rasters where lakes and coast are masked out. 
# Here we assume that the soil variables are highly correlated and fill any
# cells that are on the border of the original raster to have the mean of all 
# surrounding cells that are not NA. 

# Define function
extend_by_mean <- function(cell_values){
  # Get value of cell in the centre of the focus
  focal_cell <- cell_values[5]
  # Get values of all neighbouring cells
  other_cells <- cell_values[-5]
  # Check whether cell is NA
  if(is.na(focal_cell)){
    # cat("Cell is NA:", focal_cell, "\n")
    # If yes, check whether any of the other cells are not NA
    if(sum(is.na(other_cells)) < 8){
      # If more than one cell is not NA calculate the mean
      mean_other_cells <- mean(other_cells, na.rm = T)
      # cat("Average of other cells:", mean_other_cells, "( n = ", sum(is.na(other_cells)), ")\n")
      # Return the mean as the new value for the cell
      return(mean_other_cells)
    } else {
      # cat("All other cells are NA, returning NA.\n")
      # If all other cells are NA, return NA
      return(NA)
    }
  } else {
    # cat("Cell not NA, returning cell value:", focal_cell, "\n")
    # If the cell is not NA, then return the focal cell value
    return(focal_cell)
  }
}

# Gap fill the soil variables using focal and the helper function
clay_gap_filled <- focal(clay, fun = extend_by_mean, 
                         filename = "data/predictor_data/soil_layers/Clay_gap_filled.tif")
sand_gap_filled <- focal(sand, fun = extend_by_mean, 
                         filename = "data/predictor_data/soil_layers/sand_gap_filled.tif")
soil_carbon_gap_filled <- focal(soil_carbon, fun = extend_by_mean, 
                         filename = "data/predictor_data/soil_layers/soil_carbon_gap_filled.tif")

# Reproject (not gap-filled for training and gap-filled for projections)
project(clay, dtm10m, filename = "data/predictor_data/soil_layers/Clay_utm32_10m.tif",
        method = "near")
project(sand, dtm10m, filename = "data/predictor_data/soil_layers/Sand_utm32_10m.tif",
        method = "near")
project(soil_carbon, dtm10m, filename = "data/predictor_data/soil_layers/Soc_utm32_10m.tif",
        method = "near")
project(clay_gap_filled, dtm10m, filename = "data/predictor_data/soil_layers/gap_filled/Clay_utm32_10m.tif",
        method = "near")
project(sand_gap_filled, dtm10m, filename = "data/predictor_data/soil_layers/gap_filled/Sand_utm32_10m.tif",
        method = "near")
project(soil_carbon_gap_filled, dtm10m, filename = "data/predictor_data/soil_layers/gap_filled/Soc_utm32_10m.tif",
        method = "near")

# Set CRS for ground water raster and reproject also 
crs(ns_groundwater_summer) <- "EPSG:25832"
project(ns_groundwater_summer, dtm10m, 
        filename ="data/predictor_data/terraennaert_grundvand_10m/ns_groundwater_summer_utm32_10m.tif",
        method = "bilinear")
