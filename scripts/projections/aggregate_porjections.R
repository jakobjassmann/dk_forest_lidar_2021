# Script to generate "forest" aggregates of the 10 m data for the projections.

# Premise: The 10 m resolution data is not practical for practitioners. They 
# require larger units they can visit and inspect for value. To help with that,
# we have decided to also produce a 100 m aggregate of the projections. 
# Here we will assign a the forest quality as high quality if the percentage
# of high quality pixels in the sub-pixels exceeds the expected value of high
# quality pixels in the training polygons. 

# Dependencies
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(parallel)
library(pbapply)
library(exactextractr)
library(ggplot2)
library(cowplot)

## 1) Establish distribution of high quality pixels in high quality training polys

# Load high quality forest polygons
load("data/training_data/polygon_geometries/high_quality_polys.Rda")

## Extract proportions of high quality forest pixels

# ... for the ranger biowide projections
ranger_biowide_projections <- rast("data/projections/ranger_biowide/forest_quality_ranger_biowide.vrt")
system.time(ranger_biowide_projections_prop_high <- exact_extract(ranger_biowide_projections, high_quality))
ranger_biowide_projections_prop_high <- ranger_biowide_projections_prop_high %>%
  pblapply(function(x){
    pixel_vals <- x[,1]
    pixel_vals <- na.omit(pixel_vals)
    prop_high <- sum(pixel_vals == 1) / length(pixel_vals)  
  }) %>%
  unlist() 
ranger_biowide_projections_prop_high %>% hist() 

# ... for the ranger derek projections
ranger_derek_projections <- rast("data/projections/ranger_derek/forest_quality_ranger_derek.vrt")
system.time(ranger_derek_projections_prop_high <- exact_extract(ranger_derek_projections, high_quality))
ranger_derek_projections_prop_high <- ranger_derek_projections_prop_high %>%
  pblapply(function(x){
    pixel_vals <- x[,1]
    pixel_vals <- na.omit(pixel_vals)
    prop_high <- sum(pixel_vals == 1) / length(pixel_vals)  
  }) %>%
  unlist() 
ranger_derek_projections_prop_high %>% hist() 

# ... for the gbm biowide projections
gbm_biowide_projections <- rast("data/projections/gbm_biowide/forest_quality_gbm_biowide.vrt")
  gbm_biowide_projections_prop_high <- exact_extract(gbm_biowide_projections, high_quality)
  gbm_biowide_projections_prop_high <- gbm_biowide_projections_prop_high %>%
    pblapply(function(x){
      pixel_vals <- x[,1]
      pixel_vals <- na.omit(pixel_vals)
      prop_high <- sum(pixel_vals == 1) / length(pixel_vals)   
    }) %>%
    unlist()
  
# ... for the gbm derek projections
gbm_derek_projections <- rast("data/projections/gbm_derek/forest_quality_gbm_derek.vrt")
gbm_derek_projections_prop_high <- exact_extract(gbm_derek_projections, high_quality)
gbm_derek_projections_prop_high <- gbm_derek_projections_prop_high %>%
  pblapply(function(x){
    pixel_vals <- x[,1]
    pixel_vals <- na.omit(pixel_vals)
    prop_high <- sum(pixel_vals == 1) / length(pixel_vals)   
  }) %>%
  unlist()
  
# Identify threshold of min proportion of pixels high quality so that 80% of 
# All high quality forests are included. 
ranger_biowide_thresh <- quantile(ranger_biowide_projections_prop_high, 0.20, na.rm = T)
ranger_derek_thresh <- quantile(ranger_derek_projections_prop_high, 0.20, na.rm = T)
gbm_biowide_thresh <- quantile(gbm_biowide_projections_prop_high, 0.20, na.rm = T)
gbm_derek_thresh <- quantile(gbm_derek_projections_prop_high, 0.20, na.rm = T)


# Graph distributions of proportions
plot_grid(ggplot() +
            geom_histogram(aes(x = gbm_biowide_projections_prop_high)) +
            labs(x = "Proportion high quality pixels", y = "Count",
                 title = "GBM BIOWIDE projections") +
            scale_y_continuous(limits = c(0,3000)) +
            geom_vline(aes(xintercept = gbm_biowide_thresh), color = "red") +
            annotate("text", x = gbm_biowide_thresh, y = 3000, 
                      label = paste0(" 80% of forests\n above: ", round(gbm_biowide_thresh, 2)),
                      color = "red", hjust = 0, vjust = 1) +
            theme_cowplot(),
          ggplot() + geom_histogram(aes(x = gbm_derek_projections_prop_high)) +
            labs(x = "Proportion high quality pixels", y = "Count",
                 title = "GBM Derek projections") +
            scale_y_continuous(limits = c(0,3000)) +
            geom_vline(aes(xintercept = gbm_derek_thresh), color = "red") +
            annotate("text", x = gbm_derek_thresh, y = 3000, 
                     label = paste0(" 80% of forests\n above: ", round(gbm_derek_thresh, 2)),
                     color = "red", hjust = 0, vjust = 1) +
            theme_cowplot(),
          ggplot() +
            geom_histogram(aes(x = ranger_biowide_projections_prop_high)) +
            labs(x = "Proportion high quality pixels",  y = "Count",
                 title = "Ranger BIOWIDE projections") +
            scale_y_continuous(limits = c(0,3000)) +            
            geom_vline(aes(xintercept = ranger_biowide_thresh), color = "red") +
            annotate("text", x = ranger_biowide_thresh, y = 3000, 
                      label = paste0(" 80% of forests\n above: ", round(ranger_biowide_thresh, 2)),
                      color = "red", hjust = 0, vjust = 1) +
            theme_cowplot(),
          ggplot() +
            geom_histogram(aes(x = ranger_derek_projections_prop_high)) +
            labs(x = "Proportion high quality pixels",  y = "Count",
                 title = "Ranger DEREK projections") +
            scale_y_continuous(limits = c(0,3000)) +            
            geom_vline(aes(xintercept = ranger_derek_thresh), color = "red") +
            annotate("text", x = ranger_derek_thresh, y = 3000, 
                     label = paste0(" 80% of forests\n above: ", round(ranger_derek_thresh, 2)),
                     color = "red", hjust = 0, vjust = 1) +
            theme_cowplot(),
          ncol = 2,
          nrow = 2)

## Generate focal aggregates using the threshold

# Define aggregation functions. These function take the threshold and return
# 1 for high quality and 2 for low quality forests if the threshold is exceeded.

# Original written as one function where threshold was an argument, but these
# don't seeem to work with terra::aggregate 

forest_qual_by_threshold <- function(...){
  # Parse arguments
  arguments <- list(...)
  values <- arguments[[1]]
  if("threshold" %in% names(arguments)){
    threshold <- arguments$threshold
  } else {
    threshold <- 0.75
  } 
  
  # remove NAs 
  values <- na.omit(values)
  
  # check whether all values are NA, if so return NA
  if(length(values) == 0) return(NA)
  
  # calc proportion
  proportion <- sum((values == 1)) / length(values)
  
  # return 1 (high quality) if threshold is exceeded
  if(proportion >= threshold) return(1)
  
  # Otherwise return 2 (low quality)
  return(2)
}

# Test the function - first line should be 1, second line 2, third 2, last NA
forest_qual_by_threshold(c(rep(1, ceiling(100 * ranger_biowide_thresh)), 
                           rep(2, floor(100 * (1 - ranger_biowide_thresh)))),
                           threshold = ranger_biowide_thresh)
forest_qual_by_threshold(c(rep(1, ceiling(100 * ranger_biowide_thresh)), 
                           rep(2, 10 + floor(100 * (1 - ranger_biowide_thresh)))),
                           ranger_biowide_thresh)
forest_qual_by_threshold(c(rep(1, ceiling(100 * ranger_biowide_thresh)), 
                           rep(2, 10 + floor(100 * (1 - ranger_biowide_thresh))), 
                           rep(NA, 10)),
                         ranger_biowide_thresh)
forest_qual_by_threshold(c(NA, NA, NA),
                         threshold = ranger_biowide_thresh)
# Test the function - first line should be 1, second line 2, third 2, last NA
forest_qual_by_threshold(c(rep(1, ceiling(100 * ranger_derek_thresh)), 
                           rep(2, floor(100 * (1 - ranger_derek_thresh)))),
                         threshold = ranger_derek_thresh)
forest_qual_by_threshold(c(rep(1, ceiling(100 * ranger_derek_thresh)), 
                           rep(2, 10 + floor(100 * (1 - ranger_derek_thresh)))),
                         ranger_derek_thresh)
forest_qual_by_threshold(c(rep(1, ceiling(100 * ranger_derek_thresh)), 
                           rep(2, 10 + floor(100 * (1 - ranger_derek_thresh))), 
                           rep(NA, 10)),
                         ranger_derek_thresh)
forest_qual_by_threshold(c(NA, NA, NA),
                         threshold = ranger_derek_thresh)
# Test the function - first line should be 1, second line 2, third 2, last NA
forest_qual_by_threshold(c(rep(1, ceiling(100 * gbm_biowide_thresh)), 
                           rep(2, floor(100 * (1 - gbm_biowide_thresh)))),
                         gbm_biowide_thresh)
forest_qual_by_threshold(c(rep(1, ceiling(100 * gbm_biowide_thresh)), 
                           rep(2, 10 + floor(100 * (1 - gbm_biowide_thresh)))),
                         gbm_biowide_thresh)
forest_qual_by_threshold(c(rep(1, ceiling(100 * gbm_biowide_thresh)), 
                           rep(2, 10 + floor(100 * (1 - gbm_biowide_thresh))), 
                           rep(NA, 10)),
                         gbm_biowide_thresh)
forest_qual_by_threshold(c(NA, NA, NA),
                         threshold = gbm_biowide_thresh)
# Test the function - first line should be 1, second line 2, third 2, last NA
forest_qual_by_threshold(c(rep(1, ceiling(100 * gbm_derek_thresh)), 
                           rep(2, floor(100 * (1 - gbm_derek_thresh)))),
                         gbm_derek_thresh)
forest_qual_by_threshold(c(rep(1, ceiling(100 * gbm_derek_thresh)), 
                           rep(2, 10 + floor(100 * (1 - gbm_derek_thresh)))),
                         gbm_derek_thresh)
forest_qual_by_threshold(c(rep(1, ceiling(100 * gbm_derek_thresh)), 
                           rep(2, 10 + floor(100 * (1 - gbm_derek_thresh))), 
                           rep(NA, 10)),
                         gbm_derek_thresh)
forest_qual_by_threshold(c(NA, NA, NA),
                         threshold = gbm_derek_thresh)
# Nice one!

# Check extend is a cleanly divided by 100
(ext(ranger_biowide_projections)[2] - ext(ranger_biowide_projections)[1]) %% 100
(ext(ranger_biowide_projections)[4] - ext(ranger_biowide_projections)[3]) %% 100
(ext(ranger_derek_projections)[2] - ext(ranger_derek_projections)[1]) %% 100
(ext(ranger_derek_projections)[4] - ext(ranger_derek_projections)[3]) %% 100
(ext(gbm_biowide_projections)[2] - ext(gbm_biowide_projections)[1]) %% 100
(ext(gbm_biowide_projections)[4] - ext(gbm_biowide_projections)[3]) %% 100
(ext(gbm_derek_projections)[2] - ext(gbm_derek_projections)[1]) %% 100
(ext(gbm_derek_projections)[4] - ext(gbm_derek_projections)[3]) %% 100

# aggregate rasters to 100 m x 100 m using the threshold function:
ranger_biowide_projections_100m <- aggregate(ranger_biowide_projections, fact = 10,
          fun = forest_qual_by_threshold,
          threshold = ranger_biowide_thresh,
          filename = "data/projections/ranger_biowide/forest_quality_ranger_biowide_100m.tif",
          cores = 46, 
          overwrite = T)
ranger_derek_projections_100m <- aggregate(ranger_derek_projections, fact = 10,
                                             fun = forest_qual_by_threshold,
                                             threshold = ranger_derek_thresh,
                                             filename = "data/projections/ranger_derek/forest_quality_ranger_derek_100m.tif",
                                             cores = 46, 
                                             overwrite = T)
gbm_biowide_projections_100m <- aggregate(gbm_biowide_projections, fact = 10,
          fun = forest_qual_by_threshold,
          threshold = gbm_biowide_thresh,
          filename = "data/projections/gbm_biowide/forest_quality_gbm_biowide_100m.tif",
          cores = 46, 
          overwrite = T)
gbm_derek_projections_100m <- aggregate(gbm_derek_projections, fact = 10,
                                          fun = forest_qual_by_threshold,
                                          threshold = gbm_derek_thresh,
                                          filename = "data/projections/gbm_derek/forest_quality_gbm_derek_100m.tif",
                                          cores = 46, 
                                          overwrite = T)

# Calculate number and proportion of high quality 100 m cells
sum(values(ranger_biowide_projections_100m) == 1, na.rm = T)
sum(values(ranger_biowide_projections_100m) == 1, na.rm = T) / sum(!is.na(values(gbm_biowide_projections_100m)))
sum(values(ranger_derek_projections_100m) == 1, na.rm = T)
sum(values(ranger_derek_projections_100m) == 1, na.rm = T) / sum(!is.na(values(gbm_derek_projections_100m)))
sum(values(gbm_biowide_projections_100m) == 1, na.rm = T)
sum(values(gbm_biowide_projections_100m) == 1, na.rm = T) / sum(!is.na(values(gbm_biowide_projections_100m)))
sum(values(gbm_derek_projections_100m) == 1, na.rm = T)
sum(values(gbm_derek_projections_100m) == 1, na.rm = T) / sum(!is.na(values(gbm_derek_projections_100m)))

# Resample to 10 m and apply forest mask
ranger_biowide_projections_100m_downsampled <- resample(ranger_biowide_projections_100m, 
                                                ranger_biowide_projections,
                                                method = "near",
                                                filename = "data/projections/ranger_biowide/forest_quality_ranger_biowide_100m_downsampled.tif")
ranger_derek_projections_100m_downsampled <- resample(ranger_derek_projections_100m, 
                                                        ranger_derek_projections,
                                                        method = "near",
                                                        filename = "data/projections/ranger_derek/forest_quality_ranger_derek_100m_downsampled.tif")
gbm_derek_projections_100m_downsampled <- resample(gbm_derek_projections_100m, 
                                             gbm_derek_projections,
                                             method = "near",
                                             filename = "data/projections/gbm_derek/forest_quality_gbm_derek_100m_downsampled.tif")
gbm_derek_projections_100m_downsampled <- resample(gbm_derek_projections_100m, 
                                                     gbm_derek_projections,
                                                     method = "near",
                                                     filename = "data/projections/gbm_derek/forest_quality_gbm_derek_100m_downsampled.tif")

# # Apply the forest mask 
# forest_mask <- rast("data/basemap_forests/forest_mask.tif")
# forest_mask <- crop(forest_mask, ranger_projections_100m_downsampled)
# ranger_projections_100m_downsampled <- mask(ranger_projections_100m_downsampled, forest_mask)
# gbm_projections_100m_downsampled <- mask(gbm_projections_100m_downsampled, forest_mask)
# 
# writeRaster(ranger_projections_100m_downsampled,
#             "data/projections/ranger_biowide/forest_quality_ranger_biowide_100m_downsampled.tif",
#             overwrite = T)
# writeRaster(gbm_projections_100m_downsampled,
#             "data/projections/gbm_biowide/forest_quality_gbm_biowide_100m_downsampled.tif",
#             overwrite = T)

# Generate cloud optimised GeoTiffs
gdal_utils("warp",
           source = "data/projections/ranger_biowide/forest_quality_ranger_biowide_100m_downsampled.tif",
           destination = "data/projections/ranger_biowide/forest_quality_ranger_biowide_100m_downsampled_cog_epsg3857.tif",
           options = c(
             "-of", "COG",
             "-co", "RESAMPLING=NEAREST",
             "-co", "TILING_SCHEME=GoogleMapsCompatible",
             "-co", "COMPRESS=DEFLATE",
             "-co", "NUM_THREADS=46"
           ))
gdal_utils("warp",
           source = "data/projections/ranger_derek/forest_quality_ranger_derek_100m_downsampled.tif",
           destination = "data/projections/ranger_derek/forest_quality_ranger_derek_100m_downsampled_cog_epsg3857.tif",
           options = c(
             "-of", "COG",
             "-co", "RESAMPLING=NEAREST",
             "-co", "TILING_SCHEME=GoogleMapsCompatible",
             "-co", "COMPRESS=DEFLATE",
             "-co", "NUM_THREADS=46"
           ))
gdal_utils("warp",
           source = "data/projections/gbm_biowide/forest_quality_gbm_biowide_100m_downsampled.tif",
           destination = "data/projections/gbm_biowide/forest_quality_gbm_biowide_100m_downsampled_cog_epsg3857.tif",
           options = c(
             "-of", "COG",
             "-co", "RESAMPLING=NEAREST",
             "-co", "TILING_SCHEME=GoogleMapsCompatible",
             "-co", "COMPRESS=DEFLATE",
             "-co", "NUM_THREADS=46"
           ))
gdal_utils("warp",
           source = "data/projections/gbm_derek/forest_quality_gbm_derek_100m_downsampled.tif",
           destination = "data/projections/gbm_derek/forest_quality_gbm_derek_100m_downsampled_cog_epsg3857.tif",
           options = c(
             "-of", "COG",
             "-co", "RESAMPLING=NEAREST",
             "-co", "TILING_SCHEME=GoogleMapsCompatible",
             "-co", "COMPRESS=DEFLATE",
             "-co", "NUM_THREADS=46"
           ))

# End of File
