# Figure 2 for the DK Forest LiDAR manuscript
# Jakob J. Assmann j.assmann@bio.au.dk 16 March 2022

# Dependencies
library(tidyverse)
library(ggplot)
library(terra)
library(patchwork)
library(rasterVis)
library(rnaturalearth)
library(cowplot)

# define patchwork plot layout
plot_layout <- "
AAAB
AAAC
AAAD
EFGH
"

high_quality_col <- "#38AAB5"
low_quality_col <- "#A64304"

# Generate map of high quality forests (panel A)

# Forest projections
forest_quality <- rast("data/projections/ranger_biowide/forest_quality_ranger_biowide_100m_downsampled.tif")
# DK boundaries
denmark <- read_sf("data/stratification/biowide_georegions/DK/DK.shp") %>%
  st_transform(crs(forest_quality))

gplot(forest_quality, maxpixels = 100000) +
  geom_sf(data = denmark, inherit.aes = F, fill = NA, alpha = 0.5, size = 0.01) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = high_quality_col, high = low_quality_col,
                      na.value = NA)+
  annotate("text", 
           x = ext(forest_quality)[1] + 
             0.75 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           y = ext(forest_quality)[3] +
             0.9325 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           label = "High Quality", 
           colour = "black",
           size = 14 * 0.35,
           hjust = 0,
           vjust = 0.5) +
  annotate("rect", 
           xmin = ext(forest_quality)[1] + 
             0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           xmax = 20000 + ext(forest_quality)[1] + 
             0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           ymin = ext(forest_quality)[3] +
             0.9 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           ymax = 20000 + ext(forest_quality)[3] +
             0.9 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           color = "black",
           fill = high_quality_col) +
  annotate("text", 
           x = ext(forest_quality)[1] + 
             0.75 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           y = ext(forest_quality)[3] +
             0.8325 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           label = "Low Quality", 
           colour = "black",
           size = 14 * 0.35,
           hjust = 0,
           vjust = 0.5) +
  annotate("rect", 
           xmin = ext(forest_quality)[1] + 
             0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           xmax = 20000 + ext(forest_quality)[1] + 
             0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           ymin = ext(forest_quality)[3] +
             0.8 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           ymax = 20000 + ext(forest_quality)[3] +
             0.8 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           color = "black",
           fill = low_quality_col) +
  labs(title = "Predicted high quality forest: 2332 km²") +
  theme_map() +
  theme(legend.position = "none"
        #, panel.border = element_rect(colour = "black", fill = NA)
        )

## Zoom-in panels

# Set areas of interest
husby_klit <- st_bbox(c(xmin = 8.12327, ymin = 56.27823, 
                        xmax = 8.20458, ymax = 56.31252),
                      crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(crs(forest_quality)) %>%
  st_bbox()
mols_bjerge <- st_bbox(c(xmin = 10.49393, ymin = 56.20079, 
                         xmax = 1175515, ymax = 56.22565),
                       crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(crs(forest_quality)) %>%
  st_bbox()
frederiksdal <- st_bbox(c(xmin = 12.41060, ymin = 55.76796, 
                          xmax = 12.47054, ymax = 55.78682),
                        crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(crs(forest_quality)) %>%
  st_bbox()
bornholm <- st_bbox(c(xmin = 14.84336, ymin = 55.07429, 
                      xmax = 14.99530, ymax = 55.15600),
                    crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(crs(forest_quality)) %>%
  st_bbox()

# Adjust bounding boxes to fit same ration as big map using a helper function
adjust_bb <- function(bbox){
  ratio <-  (ext(forest_quality)[4] - ext(forest_quality)[3]) /
    (ext(forest_quality)[2] - ext(forest_quality)[1])
  # adjust bbox based on x ration 
  height <- bbox["ymax"] - bbox["ymin"]
  target_height <- (bbox["xmax"] - bbox["xmin"]) * ratio
  bbox["ymin"] <- bbox["ymin"] + 0.5 * height - 0.5 * target_height
  bbox["ymax"] <- bbox["ymin"] + target_height
  return(bbox)
}

# Check that this works as intended
plot(st_as_sfc(bornholm))
plot(st_as_sfc(adjust_bb(bornholm)), col = "red", add = T)

# Update bboxes
husby_klit <- adjust_bb(husby_klit)
mols_bjerge <- adjust_bb(mols_bjerge)
frederiksdal <- adjust_bb(frederiksdal)
bornholm <- adjust_bb(bornholm)

# helper function to generate sub panels
generate_closeup <- function(bbox){
  forest_crop <- crop(forest_quality, vect(st_as_sfc(bbox)))
  gplot(forest_crop, maxpixels = 100000) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = high_quality_col, high = low_quality_col,
                        na.value = NA) +
    coord_equal() +
    theme_map() 
}

husby_klit_forest_qual <- generate_closeup(husby_klit)
mols_bjerge_forest_qual <- generate_closeup(mols_bjerge)
frederiksdal_forest_qual <- generate_closeup(frederiksdal)
bornholm_forest_qual <- generate_closeup(bornholm)

ortho_files <- list.files("O:/Nat_Ecoinformatics/B_Read/LegacyData/Denmark/Orthophotos/SOF2016/Final_ecw/", "ECW$", full.names = T)

# helper function to generate orthophoto
get_ortho <- function(bbox){
  tile_footprints <- read_sf("F:/JakobAssmann/EcoDes-DK15_v1.1.0/tile_footprints/tile_footprints.shp")
  tile_ids <- st_intersection(tile_footprints, st_as_sfc(bbox)) %>% pull(tile_id)
  sapply(tile_ids, function(x) ortho_files[grepl(x, ortho_files)]) %>% map(rast) %>%
    reduce(merge)
}
library(stars)

read_stars(ortho_files[2334])

