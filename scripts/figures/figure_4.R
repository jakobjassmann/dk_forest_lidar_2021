# DK Forest LiDAR - Figure 4 - Disturbance figure
# Jakob J Assmann j.assmann@bio.au.dk 25 March 2022

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
disturbance_col <- "#7D3E8C"

# Forest projections
forest_quality <- rast("data/projections/ranger_biowide/forest_quality_ranger_biowide.vrt")

# DK boundaries
denmark <- read_sf("data/stratification/biowide_georegions/DK/DK.shp") %>%
  st_transform(crs(forest_quality))

# Forest disturbance
disturbance_since_2015 <- rast("data/forest_change_cs/disturbance_since_2015.tif")

# main_panel <- gplot(disturbance_since_2015, maxpixels = 500000) +
#   geom_sf(data = denmark, 
#           inherit.aes = F, 
#           fill = "#FAFAFA", 
#           colour = "#919191",
#           size = 1) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient(low = NA, high = disturbance_col,
#                       na.value = NA) +
#   annotate("text", 
#            x = ext(forest_quality)[1] + 
#              0.75 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            y = ext(forest_quality)[3] +
#              0.9325 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            label = "High Quality", 
#            colour = "black",
#            size = 14 * 0.35,
#            hjust = 0,
#            vjust = 0.5) +
#   annotate("rect", 
#            xmin = ext(forest_quality)[1] + 
#              0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            xmax = 20000 + ext(forest_quality)[1] + 
#              0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            ymin = ext(forest_quality)[3] +
#              0.9 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            ymax = 20000 + ext(forest_quality)[3] +
#              0.9 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            color = "black",
#            fill = high_quality_col) +
#   annotate("text", 
#            x = ext(forest_quality)[1] + 
#              0.75 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            y = ext(forest_quality)[3] +
#              0.8325 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            label = "Low Quality", 
#            colour = "black",
#            size = 14 * 0.35,
#            hjust = 0,
#            vjust = 0.5) +
#   annotate("rect", 
#            xmin = ext(forest_quality)[1] + 
#              0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            xmax = 20000 + ext(forest_quality)[1] + 
#              0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            ymin = ext(forest_quality)[3] +
#              0.8 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            ymax = 20000 + ext(forest_quality)[3] +
#              0.8 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            color = "black",
#            fill = low_quality_col) +  
#   annotate("text", 
#            x = ext(forest_quality)[1] + 
#              0.75 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            y = ext(forest_quality)[3] +
#              0.7325 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            label = "Disturbed Since 2015", 
#            colour = "black",
#            size = 14 * 0.35,
#            hjust = 0,
#            vjust = 0.5) +
#   annotate("rect", 
#            xmin = ext(forest_quality)[1] + 
#              0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            xmax = 20000 + ext(forest_quality)[1] + 
#              0.69 * (ext(forest_quality)[2] - ext(forest_quality)[1]),
#            ymin = ext(forest_quality)[3] +
#              0.7 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            ymax = 20000 + ext(forest_quality)[3] +
#              0.7 * (ext(forest_quality)[4] - ext(forest_quality)[3]),
#            color = "black",
#            fill = disturbance_col) +
#   labs(title = paste0("Disturbed High Quality Forest: 19 km²\n",
#                       "Disturbed Low Quality Forest: 54 km²")) +
#   theme_map() +
#   theme(legend.position = "none",
#         plot.margin = unit(c(0.1,0,0,0), "in"),
#         #panel.border = element_rect(colour = "red", fill = NA)
#   )


mols_bjerge <- st_bbox(c(xmin = 593625, ymin = 6230433, 
                         xmax = 593625 + 1500, ymax = 6230433 + 1000),
                       crs =  crs(forest_quality)) %>%
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
# Update bbox
mols_bjerge <- adjust_bb(mols_bjerge)

# Load orthophoto file names and tile_footprints
ortho_files_2021 <- list.files("F:/JakobAssmann/orthos_2021/extractions", "[eE][cC][wW]$", full.names = T)
ortho_files_2014 <- list.files("O:/Nat_Ecoinformatics/B_Read/LegacyData/Denmark/Orthophotos/SOF2014/UTM32N", "[eE][cC][wW]$", full.names = T)
tile_footprints <- read_sf("F:/JakobAssmann/EcoDes-DK15_v1.1.0/tile_footprints/tile_footprints.shp")

# helper function to generate orthophoto
get_ortho <- function(bbox, area_name, scale_to = 0.5, orthos = ortho_files_2014){
  tile_ids <- st_intersection(tile_footprints, st_as_sfc(bbox)) %>% pull(tile_id)
  # Get list of ecw files
  if(grepl(".*2014.*", orthos[1])) {
    tile_ids <- tile_ids %>% 
    sapply(function(tile_id){
      x <- gsub("([0-9]{4})_([0-9]{3})", "\\1", tile_id) %>% as.numeric()
      y <- gsub("([0-9]{4})_([0-9]{3})", "\\2", tile_id) %>% as.numeric()
      if(x %% 2 != 0) x <- x - 1
      if(y %% 2 != 0) y <- y - 1
      return(paste0(x, "_", y))
    }) %>%
    unique()
    }
  ecw_files <- sapply(tile_ids, function(x) orthos[grepl(x, orthos)]) 
  
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
# mols_bjerge_forest_ortho_2014 <- get_ortho(mols_bjerge, "mols_bjerge_2014",
#                                            orthos = ortho_files_2014)
# mols_bjerge_forest_ortho_2021 <- get_ortho(mols_bjerge, "mols_bjerge_2021",
#                                            orthos = ortho_files_2021)

mols_bjerge_forest_ortho_2021 <- rast("data/orthophotos/mols_bjerge_2021.tif")
mols_bjerge_forest_ortho_2014 <- rast("data/orthophotos/mols_bjerge_2014.tif")

# function to generate plot
plot_ortho_n_qual_dist <- function(ortho, ortho_name, qual = T, dist = F,
                                   dist_col = "black"){
  forest_qual_crop <- crop(forest_quality, ext(ortho))
  dist_crop <- crop(disturbance_since_2015, ext(ortho)) %>%
    as.polygons() %>%
    st_as_sf() %>%
    st_union() %>%
    vect()
  width <- ext(ortho)[2] - ext(ortho)[1]
  height <- ext(ortho)[4] - ext(ortho)[3]
  temp_file <- tempfile()
  png(temp_file, width = 772 * 2, height = 603 * 2)
  plotRGB(ortho)
  if(qual == T){
  plot(forest_qual_crop,
       col = c(high_quality_col, low_quality_col),
       alpha = 0.6,
       add = T)
  }
  if(dist == T){
    polys(dist_crop,
         col =  "#00000000",
         density = 5,
         alpha = 1,
         lwd = 10,
         #col = dist_col,
         border = dist_col)
    # polys(dist_crop,
    #       #col =  "#00000000",
    #       density = 5,
    #       alpha = 1,
    #       lwd = 2,
    #       col = dist_col,
    #       border = NA)
  }
  text(ext(ortho)[1] + width * 0.04,
       ext(ortho)[3] + height * 0.90,
       ortho_name,
       adj = 0,
       col = "white",
       cex = 10)
  rect(ext(ortho)[1] + width * 0.9 - 300,
       ext(ortho)[3] + height * 0.08,
       ext(ortho)[1] + width * 0.9,
       ext(ortho)[3] + height * 0.08 + height * 0.025,
       col = "white",
       border = "white")
  text(ext(ortho)[1] + width * 0.9 - 150,
       ext(ortho)[3] + height * 0.17,
       "300 m",
       col = "white",
       cex = 10)
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

mols_bjerge_grob_2014_qual <- plot_ortho_n_qual_dist(mols_bjerge_forest_ortho_2014,
                                                "Forest Quality", qual = T, dist = F, 
                                                dist_col = "#C575D9")

mols_bjerge_grob_2021_dist <- plot_ortho_n_qual_dist(mols_bjerge_forest_ortho_2021,
                                                     "Disturbance", qual = F, dist = T, 
                                                     dist_col = "#C575D9")

mols_bjerge_grob_2014 <- plot_ortho_n_qual_dist(mols_bjerge_forest_ortho_2014,
                                                "2014 Summer", qual = F, dist = F, 
                                                dist_col = "#C575D9")

mols_bjerge_grob_2021 <- plot_ortho_n_qual_dist(mols_bjerge_forest_ortho_2021,
                                      "2021 Spring", qual = F, dist = F,
                                      dist_col = "#C575D9")
plot_grid(mols_bjerge_grob_2014,
          mols_bjerge_grob_2014_qual,
          mols_bjerge_grob_2021,
          mols_bjerge_grob_2021_dist,
          nrow = 2,
          rel_widths = c(1,1)) %>%
  save_plot("docs/figures/figure_4.png", .,
            base_height = 8,
            base_asp = 772/603, #  (772 + 0.5 * 772) / 603,
            bg = "white")
