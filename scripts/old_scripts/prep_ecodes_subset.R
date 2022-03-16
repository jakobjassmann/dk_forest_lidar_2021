# DK-LiDAR Example script - dataset subsetting
# Jakob J. Assmann j.assmann@bio.au.dk June 2021
# Updated for Aarhus region on 9 July 2021

# Dependencies
library(sf)
library(tidyverse)
library(raster)

# Set area of interest (here Husby Klit Dune reserve)
aoi <- st_polygon(list(matrix(c(10.284891, 56.040591, 
                                10.284891,56.306278,
                                9.876908,56.306278,
                                9.876908,56.040591,
                                10.284891, 56.040591), 
                              ncol = 2,
                              byrow = T))) %>%
  st_sfc(crs = 4326)

# Read in tile footprints
tile_footprints <- read_sf("F:/JakobAssmann/EcoDes-DK_v1.1.0/tile_footprints/tile_footprints.shp")

# Determine intersecting tiles
aoi_tiles <- aoi %>%
  st_transform(., st_crs(tile_footprints)) %>%
  st_intersects(tile_footprints, ., sparse = F) %>%
  filter(tile_footprints, .)

# Get dirs for the available variables 
# (change according to EcoDes-DK15 directory, using a toy directory here)
ecodes_dir <- "F:/JakobAssmann/EcoDes-DK_v1.1.0/" 
variable_dirs <- shell(paste0("dir /b /s /a:d ", 
                              gsub("/", "\\\\", ecodes_dir)), intern = T) %>%
  gsub("\\\\", "/", .)

# Remove folders with subfolders
variable_dirs <- variable_dirs[variable_dirs %>% 
                                 map(function(x) sum(grepl(x, variable_dirs)) == 1) %>%
                                 unlist()]
variable_dirs <- variable_dirs[!grepl("tile_footprints", variable_dirs)]

# Set list of variables to exclude from subset
variables_to_exclude <- c("point_source_info.*", "point_count.*", "masks$",
                         "vegetation_proportion$", "tile_footprints")

# Subsets dirs using regex matching
sub_dirs_to_export <- variable_dirs
for(i in seq_along(variables_to_exclude)){
  sub_dirs_to_export <- sub_dirs_to_export[!grepl(variables_to_exclude[i], 
                                                  sub_dirs_to_export)]
}

# Set target dir to place output:
target_dir <- "data/ecodes_subset"
dir.create(target_dir)

# Create absolute folder paths
in_dirs <- sub_dirs_to_export
out_dirs <- paste0(target_dir, "/", gsub(ecodes_dir, "", sub_dirs_to_export))

# Define helper function to copy tiles 
copy_subset <- function(tiles, in_dir, out_dir){
  # Get variable name
  var_name <- gsub(".*/(.*)$", "\\1", in_dir)
  # Status
  cat(rep("#", 80), "\n", sep = "")
  cat("Copying", length(tiles), "tiles for variable:", var_name, "\n")
  cat("from:\t", in_dir, "\n")
  cat("to:\t", in_dir, "\n")
  # Check whether out dir exists
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
  # Copy tiles
  tiles %>% map(function(tile_id){
    file_name <- paste0(var_name, "_", tile_id, ".tif")
    in_path <- paste0(in_dir, "/", file_name)
    out_path <- paste0(out_dir, "/", file_name)
    file.copy(in_path, out_path)
    cat(".")
  })
  cat("\nDone.\n")
  return(NULL)
}

# Copy / extract subset of dataset by applying helper function
map2(in_dirs,
     out_dirs,
     function(in_dir, out_dir) copy_subset(aoi_tiles$tile_id, in_dir, out_dir))

# Generate vrt files for subset variables using sf's inbuild gdalutils 
map(out_dirs,
    function(dir_name) {
      oldwd <- getwd()
      setwd(dir_name)
      var_name <- gsub(".*/(.*)$", "\\1", dir_name)
      cat("Generating vrt for:", var_name, "\n")
      gdal_utils("buildvrt",
                 source = list.files(getwd(),".tif", full.names = T),
                 destination = paste0(var_name, ".vrt"))
      setwd(oldwd)
    })

# End of file