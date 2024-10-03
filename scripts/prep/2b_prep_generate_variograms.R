# Sample variograms for predictor rasters
# Jakob Assmann j.assmann@bio.au.dk 8 February 2022

## 1) Housekeeping ----

# Dependencies
library(raster)
library(terra)
library(gstat)
library(parallel)
library(pbapply)
library(cowplot)
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)

# set global raster options with process bar
rasterOptions(progress = "text")
terraOptions(progress=1)

# Set raster file names
raster_files <- shell("dir /b /s F:\\JakobAssmann\\EcoDes-DK15_v1.1.0\\*.vrt",
                      intern = T) %>%
  gsub("\\\\", "/", .)

# Remove all variables not relevant
raster_files <- raster_files[!grepl("date_stamp",raster_files)]
raster_files <- raster_files[!grepl("point_count",raster_files)]
raster_files <- raster_files[!grepl("point_source",raster_files)]
raster_files <- raster_files[!grepl("building",raster_files)]

# Add non-EcoDes variables
raster_files <- c(
  raster_files,
  "data/predictor_data/terraennaert_grundvand_10m/ns_groundwater_summer_utm32_10m.tif",
  "data/predictor_data/soil_layers/Clay_utm32_10m.tif",
  "data/predictor_data/soil_layers/Sand_utm32_10m.tif",
  "data/predictor_data/soil_layers/Soc_utm32_10m.tif",
  "data/predictor_data/foliage_height_diversity/foliage_height_diversity.tif")

# Remove superfluous variables
raster_files <- c("mask") %>% 
  paste0("(.*", ., ".*)", collapse = "|") %>%
  grepl(., raster_files) %>%
  `!` %>%
  raster_files[.]

# Define function to calculate variograms
sample_variograms <- function(predictor_raster_file) { 
  
  # Load raster
  predictor_raster <- raster(predictor_raster_file)
  
  # # Check whether raster is part of the EcoDes-DK dataset
  # # If so apply conversion factor
  # if(grepl("EcoDes", predictor_raster_file)){
  #   conversion_factors <- read.csv("F:/JakobAssmann/EcoDes-DK15_v1.1.0/conversion_factors.csv")
  #   conversion_factor <- conversion_factors$conv_fac[sapply(conversion_factors[,1], grepl, x = predictor_raster_file)]
  #   predictor_raster <- predictor_raster / conversion_factor
  # }
  
  # Convert raster to spdf
  predictor_spdf <- as(predictor_raster, "SpatialPixelsDataFrame" ) 
  
  # Square out spdf (needed due to a bug in gstat)
  predictor_spdf@grid@cellsize[1] <- as.numeric(formatC(predictor_spdf@grid@cellsize[1], 
                                                        format = "d"))
  predictor_spdf@grid@cellsize[2] <- as.numeric(formatC(predictor_spdf@grid@cellsize[2], 
                                                        format = "d"))
  # Autogenerate variogram forumla 
  vario_formula <- as.formula(paste0(names(predictor_raster),
                                     " ~ 1"))
  # Sample variograms (this can take ages)
  vario_1km <- variogram(vario_formula, 
                         predictor_spdf[sample(nrow(predictor_spdf) / 1000),],
                         width = 10,
                         cutoff = 1000,
                         verbose = T)
  vario_10km <- variogram(vario_formula, 
                          predictor_spdf[sample(nrow(predictor_spdf) / 1000),],
                          width = 100,
                          cutoff = 10000,
                          verbose = T)
  
  # Change id colum
  vario_1km$id <- names(predictor_raster)
  vario_10km$id <- names(predictor_raster)
  
  # save variograms
  save(vario_1km, file = paste0("data/variograms/", names(predictor_raster), "_1km.Rda"))
  save(vario_10km, file = paste0("data/variograms/", names(predictor_raster), "_10km.Rda"))
  
  # clean memory
  gc()
  
  # Return variograms
  return(list(vario_1km = vario_1km, vario_10km = vario_10km))
}


# Prep parallel environment
cl <- makeCluster(14)
clusterEvalQ(cl, {
  library(gstat)
  library(raster)
  })

## 2) Sampling of variograms

# Sample variograms
vario_list <- pblapply(raster_files, sample_variograms, cl = cl)
save(vario_list, file = "data/variograms/variogram_list.Rda")
#load("data/variograms/variogram_list.Rda")
stopCluster(cl)

# Assign names
names(vario_list) <- gsub(".*/(.*)\\..*", "\\1", raster_files)

# Plot Variograms
plot_variogram <- function(vario, raster_name){
  vario_plot_1km <- ggplot(
    vario[[1]], 
    aes(x = dist, y = gamma)) + 
    geom_point() +
    labs(subtitle = paste0(raster_name, 
                           "\n1 km, 10 m bins\nminimum samples per bin: ",
                           min(vario[[1]]$np)),
         x = "Distance (m)", y = "Semivariance") +
    theme_cowplot(15)
  vario_plot_10km <- ggplot(
    vario[[2]], 
    aes(x = dist, y = gamma)) + 
    geom_point() +
    labs(subtitle = paste0(raster_name, 
                           "\n10 km, 100 m bins\nminimum samples per bin: ",
                           min(vario[[2]]$np)),
         x = "Distance (m)", y = "Semivariance") +
    theme_cowplot(15)
  vario_plots <- plot_grid(vario_plot_1km,
            vario_plot_10km, 
            ncol = 2)
  save_plot(paste0("docs/variograms/", raster_name, ".png"),
            vario_plots,
            ncol = 2,
            bg = "white")
  return("OK")
}

# Plot variograms
mapply(plot_variogram, vario_list, names(vario_list))

# End of file