# Figure 2 for the DK Forest LiDAR manuscript
# Jakob J. Assmann j.assmann@bio.au.dk 16 March 2022

# Dependencies
library(tidyverse)
library(ggplot2)
library(terra)
library(sf)
library(patchwork)
library(rasterVis)
library(rnaturalearth)
library(cowplot)

# Define colours
high_quality_col <- "#3CB0AE"
low_quality_col <- "#D67D49" 

# Generate map of high quality forests (panel A)

# Forest projections
forest_quality <- rast("data/projections/ranger_biowide/forest_quality_ranger_biowide.vrt")
# DK boundaries
denmark <- read_sf("data/stratification/biowide_georegions/DK/DK.shp") %>%
  st_transform(crs(forest_quality))

main_panel <- gplot(forest_quality, maxpixels = 500000) +
  geom_sf(data = denmark, 
          inherit.aes = F, 
          fill = "#FAFAFA", 
          colour = "#919191",
          size = 0.5) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = high_quality_col, high = low_quality_col,
                      na.value = NA)+
  annotate("text", 
           x = ext(forest_quality)[1] + 
             0.85 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           y = ext(forest_quality)[3] +
             0.9325 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           label = "High Value", 
           colour = "black",
           size = 14 * 0.35,
           hjust = 0,
           vjust = 0.5) +
  annotate("rect", 
           xmin = ext(forest_quality)[1] + 
             0.79 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           xmax = 20000 + ext(forest_quality)[1] + 
             0.79 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           ymin = ext(forest_quality)[3] +
             0.9 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           ymax = 20000 + ext(forest_quality)[3] +
             0.9 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           color = "black",
           fill = high_quality_col) +
  annotate("text", 
           x = ext(forest_quality)[1] + 
             0.85 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           y = ext(forest_quality)[3] +
             0.8325 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           label = "Low Value", 
           colour = "black",
           size = 14 * 0.35,
           hjust = 0,
           vjust = 0.5) +
  annotate("rect", 
           xmin = ext(forest_quality)[1] + 
             0.79 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           xmax = 20000 + ext(forest_quality)[1] + 
             0.79 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           ymin = ext(forest_quality)[3] +
             0.8 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           ymax = 20000 + ext(forest_quality)[3] +
             0.8 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
           color = "black",
           fill = low_quality_col) +
  annotate("rect",
           xmin = ext(forest_quality)[1] + 
             0.675 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           xmax = 100000 + ext(forest_quality)[1] + 
             0.675 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           ymin = st_bbox(denmark)["ymin"] + 20000,
           ymax = st_bbox(denmark)["ymin"] + 20000 + 5000,
           fill = "black") +
  annotate("text", 
           x = 50000 + ext(forest_quality)[1] + 
             0.675 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
           y = st_bbox(denmark)["ymin"] + 20000 + 15000,
           label = "100 km", 
           colour = "black",
           size = 14 * 0.35,
           hjust = 0.5,
           vjust = 0.5) +
  labs(title = "Estimated High Conservation Value Forest: 1420 km² (141958 ha)") +
  theme_map() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0,0,0), "in"),
        #panel.border = element_rect(colour = "red", fill = NA)
        )

## Zoom-in panels

# Set areas of interest
husby_klit <- st_bbox(c(xmin = 445019, ymin = 6237933, 
                        xmax = 445019 + 3500, ymax = 6237933 + 3000),
                      crs = crs(forest_quality)) %>%
  st_as_sfc() %>%
  st_bbox()
mols_bjerge <- st_bbox(c(xmin = 592625, ymin = 6229433, 
                         xmax = 592625 + 3500, ymax = 6229433 + 3000),
                       crs =  crs(forest_quality)) %>%
  st_as_sfc() %>%
  st_bbox()
frederiksdal <- st_bbox(c(xmin = 713834, ymin = 6185156, 
                          xmax = 713834 + 3500, ymax = 6185156 + 3000),
                        crs =  crs(forest_quality)) %>%
  st_as_sfc() %>%
  st_bbox()
bornholm <- st_bbox(c(xmin = 877137, ymin = 6122035, 
                      xmax = 877137 + 3500, ymax = 6122035 + 3000),
                    crs = crs(forest_quality)) %>%
  st_as_sfc() %>%
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

# Load orthophoto file names and tile_footprints
# ortho_files <- list.files("O:/Nat_Ecoinformatics/B_Read/LegacyData/Denmark/Orthophotos/SOF2014/UTM32N", "[eE][cC][wW]$", full.names = T)
# tile_footprints <- read_sf("F:/JakobAssmann/EcoDes-DK15_v1.1.0/tile_footprints/tile_footprints.shp")

# helper function to generate orthophoto
get_ortho <- function(bbox, area_name, scale_to = 0.5){
  tile_ids <- st_intersection(tile_footprints, st_as_sfc(bbox)) %>% pull(tile_id)
  # Get list of ecw files
  tile_ids <- tile_ids %>% 
    sapply(function(tile_id){
      x <- gsub("([0-9]{4})_([0-9]{3})", "\\1", tile_id) %>% as.numeric()
      y <- gsub("([0-9]{4})_([0-9]{3})", "\\2", tile_id) %>% as.numeric()
      if(x %% 2 != 0) x <- x - 1
      if(y %% 2 != 0) y <- y - 1
      return(paste0(x, "_", y))
    }) %>%
    unique()
  ecw_files <- sapply(tile_ids, function(x) ortho_files[grepl(x, ortho_files)]) 
  
  # Set scale term if needed
  if(scale_to == FALSE){
    scale_to <- ""
  } else{
    scale_to <- paste0("-tr ", scale_to, " ", scale_to, " ")
  }
  
  # Convert to GeoTiff
  ecw_files %>%
    map(function(x){
      if(!file.exists(paste0("data/orthophotos/", gsub(".*/(.*)\\.[eE][cC][wW]", "\\1", x), ".tif")))
      shell(paste0("C:/OSGeo4W/OSGeo4W.bat gdalwarp ",
                   scale_to,
                   "-co COMPRESS=DEFLATE ",
                   x,
                 " data/orthophotos/", gsub(".*/(.*)\\.[eE][cC][wW]", "\\1", x), ".tif"))
    })
  geo_tiffs <- paste0("data/orthophotos/", gsub(".*/(.*)\\.[eE][cC][wW]", "\\1", ecw_files), ".tif")
  # If needed merge write out and remove tiles
  if(length(geo_tiffs) > 1) {
    geo_tiffs %>% 
    map(rast) %>% 
    reduce(merge) %>%
    crop(vect(st_as_sfc(bbox))) %>%
    writeRaster(., paste0("data/orthophotos/", area_name, ".tif"), gdal = c("COMPRESS=DEFLATE"))
    file.remove(geo_tiffs)
  } else {
    rast(geo_tiffs) %>% 
      crop(vect(st_as_sfc(bbox))) %>%
      writeRaster(., paste0("data/orthophotos/", area_name, ".tif"), gdal = c("COMPRESS=DEFLATE"))
    file.remove(geo_tiffs)
  }
  return(rast(paste0("data/orthophotos/", area_name, ".tif")))
}

# Generate orthos
# husby_klit_forest_ortho <- get_ortho(husby_klit, "husby_klit")
# mols_bjerge_forest_ortho <- get_ortho(mols_bjerge, "mols_bjerge")
# frederiksdal_forest_ortho <- get_ortho(frederiksdal, "frederiksdal")
# bornholm_forest_ortho <- get_ortho(bornholm, "bornholm")

# husby_klit_forest_ortho <- rast("data/orthophotos/husby_klit.tif")
# mols_bjerge_forest_ortho <- rast("data/orthophotos/mols_bjerge.tif")
frederiksdal_forest_ortho <- rast("data/orthophotos/frederiksdal.tif")
# bornholm_forest_ortho <- rast("data/orthophotos/bornholm.tif")

plotRGB(frederiksdal_forest_ortho,
        scale = 255)
# function to generate plot
plot_ortho_n_qual <- function(ortho, ortho_name, qual = F){
  forest_qual_crop <- crop(forest_quality, ext(ortho))
  width <- ext(ortho)[2] - ext(ortho)[1]
  height <- ext(ortho)[4] - ext(ortho)[3]
  temp_file <- tempfile()
  png(temp_file, width = 772 * 2, height = 603 * 2)
  plotRGB(ortho)
  if(qual == T){
    plot(forest_qual_crop,
         col = c(high_quality_col, low_quality_col),
         alpha = 0.8,
         add = T)
  } 
  text(ext(ortho)[1] + width * 0.04,
       ext(ortho)[3] + height * 0.88,
       ortho_name,
       adj = 0,
       col = "white",
       cex = 10)
  rect(ext(ortho)[1] + width * 0.9 - 1000,
       ext(ortho)[3] + height * 0.1,
       ext(ortho)[1] + width * 0.9,
       ext(ortho)[3] + height * 0.1 + height * 0.025,
       col = "white",
       border = "white")
  text(ext(ortho)[1] + width * 0.9 - 500,
       ext(ortho)[3] + height * 0.18,
       "1 km",
       col = "white",
       cex = 8)
  dev.off()
  gg_grob <- ggplot() +
    draw_image(temp_file) +
    geom_rect(aes(xmin = 0, xmax = 1,
                  ymin = 0.5 - (0.5 * (603) / (772)), 
                  ymax = 0.5 + (0.5 * (603) / (772))), 
              colour = "black",
              size = 0.5,
              fill = NA) +
    # labs(subtitle = ortho_name) +
    theme_map()  +
    theme(plot.margin=unit(c(0,0,0,0),"mm")
          #panel.border = element_rect(colour = "red", fill = NA)
          )
  rm(temp_file)
  return(gg_grob)
}

# husby_klit_grob <- plot_ortho_n_qual(husby_klit_forest_ortho,
#                                      "Husby Klit")
# mols_bjerge_grob <- plot_ortho_n_qual(mols_bjerge_forest_ortho,
#                                       "Mols Bjerge")
frederiksdal_grob <- plot_ortho_n_qual(frederiksdal_forest_ortho,
                                       "Frederiksdal")
frederiksdal_qual_grob <- plot_ortho_n_qual(frederiksdal_forest_ortho,
                                       "Frederiksdal", qual = T)
# bornholm_grob <- plot_ortho_n_qual(bornholm_forest_ortho,
#                                    "Bornholm")


# Compose plots

plot_grid(main_panel,
          plot_grid(frederiksdal_grob,
                    frederiksdal_qual_grob,
                    nrow = 2),
          ncol = 2,
          rel_widths = c(2,1)) %>%
  save_plot("docs/figures/figure_3.png", .,
            base_height = 6,
            base_asp = (772 + 0.5 * 772) / 603,
            bg = "white")
