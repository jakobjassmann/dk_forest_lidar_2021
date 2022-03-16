# Script to generate a foliage height diversity raster from EcoDes-DK15 
# Jakob J. Assmann j.assmann@bio.au.dk 17 February 2022

# Dependencies
library(terra)
terraOptions(progress = 2)

# Load list of EcoDes-DK15 descriptors
raster_list <- shell("dir /b /s F:\\JakobAssmann\\EcoDes-DK_v1.1.0\\*.vrt",
                     intern = T) 
raster_list <- gsub("\\\\", "/", raster_list)

# Filter out all non-proportion rasters
raster_list <- raster_list[grepl("vegetation_proportion", raster_list)]

# Stratify vegetation according to Wilson (1974): 0–1.5 m, 1.5–9 m, and >9 m
prop_below_1.5 <- raster_list[1:3]
prop_below_9 <- raster_list[4:11]
prop_above_9 <- raster_list[12:24]

# Calculate cumulative proportions for layers
prop_below_1.5 <- sum(rast(prop_below_1.5))
prop_below_9 <- sum(rast(prop_below_9))
prop_above_9 <- sum(rast(prop_above_9))

# Adjust for rounding inaccuracies
prop_below_1.5[prop_below_1.5 > 10000] <- 10000
prop_below_9[prop_below_9 > 10000] <- 10000
prop_above_9[prop_above_9 > 10000] <- 10000

# Apply correction factor to convert to actual proportion
prop_below_1.5 <- prop_below_1.5 / 10000
prop_below_9 <- prop_below_9 / 10000
prop_above_9 <- prop_above_9 / 10000

# Calculate log proportion
prop_below_1.5_log <- log(prop_below_1.5)
prop_below_9_log <- log(prop_below_9)
prop_above_9_log <- log(prop_above_9)

# Set log proportion to 0 if proportion is 0
# This is done by convention as log 0 is not defined.
# (see https://stats.stackexchange.com/questions/57069/alternative-to-shannons-entropy-when-probability-equal-to-zero)
prop_below_1.5_log[prop_below_1.5 == 0] <- 0
prop_below_9_log[prop_below_9 == 0] <- 0
prop_above_9_log[prop_above_9 == 0] <- 0

# Calculate foliage height diversity -SUM(pi x log(pi))
foliage_height_div <- -1 * ((prop_below_1.5 * prop_below_1.5_log) +
                              (prop_below_9 * prop_below_9_log) +
                              (prop_above_9 * prop_above_9_log))

# Write out raster
writeRaster(foliage_height_div,
            "data/predictor_data/foliage_height_diversity/foliage_height_diversity.tif")

# Quality control calculations
test_points <- vect(matrix(c(533766.616,6219861.380,
  529345.852,6219203.933,
  874920.332,6130082.459),
  byrow = T, ncol = 2),
  crs = crs(foliage_height_div))

# Extract qc data
qc_data <- data.frame(
  prop_below_1.5 = extract(prop_below_1.5, test_points)[,2],
  prop_below_9 = extract(prop_below_9, test_points)[,2],
  prop_above_9 = extract(prop_above_9, test_points)[,2],
  prop_below_1.5_log = extract(prop_below_1.5_log, test_points)[,2],
  prop_below_9_log = extract(prop_below_9_log, test_points)[,2],
  prop_above_9_log = extract(prop_above_9_log, test_points)[,2],
  foliage_height_div = extract(foliage_height_div, test_points)[,2]
)
qc_data$foliage_height_dif_qc <- -1 * 
  ((qc_data$prop_below_1.5 * qc_data$prop_below_1.5_log) + 
     (qc_data$prop_below_9 * qc_data$prop_below_9_log) +
     (qc_data$prop_above_9 * qc_data$prop_above_9_log))
# Looking good! :)

# EOF