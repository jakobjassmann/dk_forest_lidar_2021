# DK Forest LiDAR manuscript figure 1 - map
# Brief script to generate the map figure used in Figure 1 - workflow
# Jakob J. Assmann j.assmann@bio.au.dk 5 April 2022

# Dependencies
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(ggtext)

# Load shape files 
biowide_regions <- read_sf("data/stratification/biowide_georegions/biowide_zones.shp")

biowide_regions <- mutate(
  biowide_regions,
  x_cen = st_coordinates(st_centroid(biowide_regions))[,1],
  y_cen = st_coordinates(st_centroid(biowide_regions))[,2],
  map_span_x = st_bbox(biowide_regions)["xmax"]- 
    st_bbox(biowide_regions)["xmin"],
  map_span_y = st_bbox(biowide_regions)["ymax"]- 
    st_bbox(biowide_regions)["ymin"]) 

biowide_regions <- 
  mutate(biowide_regions,
         x = x_cen + map_span_x * c(-0.10,  # Nordjlland
                                    -0.15,  # Vestjylland
                                    0.30,  # Oestjylland
                                    0.20,  # Sjaelland
                                    -0.00,  # Bornholm
                                    -0.20), # Fune_Lolland
         
         y = y_cen + map_span_y * c( 0.10,  # Nordjlland
                                     0.10,  # Vestjylland
                                     0.30,  # Oestjylland
                                     0.10,  # Sjaelland
                                     -0.10,  # Bornholm
                                     -0.20), # Fune_Lolland
         hjust = c(1, # Nordjlland
                   1, # Vestjylland
                   0, # Oestjylland
                   0, # Sjaelland
                   0, # Bornholm
                   0),# Fune_Lolland
         vjust = c(0, # Nordjlland
                   1, # Vestjylland
                   0, # Oestjylland
                   0, # Sjaelland
                   1, # Bornholm
                   1) # Fune_Lolland
  )
biowide_regions$region[6] <- "Fune-Lolland"

# Set map extent
main_panel_xlim <- st_bbox(biowide_regions)[c(1,3)] * c(0.6, 1.15)
main_panel_ylim <- st_bbox(biowide_regions)[c(2,4)] * c(0.98, 1.0125)
main_panel_width <- main_panel_xlim[2] - main_panel_xlim[1]
main_panel_height <- main_panel_ylim[2] - main_panel_ylim[1]

# Plot map
(map_plot <- ggplot() + 
    geom_sf(aes(colour = region,
                fill = region), 
            size = 0.5,
            data = biowide_regions) +
    geom_richtext(aes(x = x, 
                      y = y, 
                      label = paste0("<span style = font-size:16pt>", region, "</span>"),
                      colour = region,
                      hjust = hjust,
                      vjust = vjust), 
                  data = biowide_regions,
                   fill = NA,
                  label.color = NA
              ) +
    scale_colour_manual(values = c("#0F403F", # Bornholm 
                                   "#3D8A88", # Fune_Lolland 
                                   "#C575D9", # Nordjlland 
                                   "#7D3E8C", # Oestjylland 
                                   "#62B5B4", # Sjaelland 
                                   "#B88AC5")) +  # Vestjylland
    scale_fill_manual(values = c("#F5FAFA", # Bornholm
                                 "#F5FAFA", # Fune_Lolland
                                 "#F9F5FA", # Nordjlland
                                 "#F9F5FA", # Oestjylland
                                 "#F5FAFA", # Sjaelland
                                 "#F9F5FA")) + # Vestjylland
    coord_sf(clip = "off",
             crs = st_crs(biowide_regions),
             xlim = main_panel_xlim,
             ylim = main_panel_ylim) +
    theme_map() +
    theme(legend.position = "none",
          #plot.background = element_rect(color = "red")
          ))
save_plot("docs/figures/map_for_figure_1.png",
          map_plot,
          base_asp = main_panel_width / main_panel_height,
          bg = "white")
  
  

