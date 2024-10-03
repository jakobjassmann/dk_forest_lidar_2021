# DK Forest LiDAR manuscript figure 1 - map
# Brief script to generate the map figure used in Figure 1 - workflow
# Jakob J. Assmann j.assmann@bio.au.dk 5 April 2022

# Dependencies
library(tidyverse)
library(sf)
library(ggplot2)
library(cowplot)
library(ggtext)
library(rnaturalearth)

# Load shape files 
biowide_regions <- read_sf("data/stratification/biowide_georegions/biowide_zones.shp")
biowide_regions$region <- factor(
  c("Northern<br>Jutland", 
  "Western<br>Jutland", 
  "Eastern<br>Jutland", 
  "Zealand", 
  "Bornholm", 
  "Funen-Lolland<br>Falster-Møn"),
  levels = 
  c("Northern<br>Jutland", 
  "Western<br>Jutland", 
  "Eastern<br>Jutland", 
  "Zealand", "Bornholm", 
  "Funen-Lolland<br>Falster-Møn")
)

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
         x = x_cen + map_span_x * c(-0.3,  # Nordjlland
                                    -0.35,  # Vestjylland
                                    0.5,  # Oestjylland
                                    0.15,  # Sjaelland
                                    -0.05,  # Bornholm
                                    -0.0), # Fune_Lolland
         
         y = y_cen + map_span_y * c( 0.05,  # Nordjlland
                                     0.10,  # Vestjylland
                                     0.30,  # Oestjylland
                                     0.05,  # Sjaelland
                                     -0.10,  # Bornholm
                                     -0.15), # Fune_Lolland
         hjust = c(0.5, # Nordjlland
                   0.5, # Vestjylland
                   0.5, # Oestjylland
                   0, # Sjaelland
                   0, # Bornholm
                   0.5),# Fune_Lolland
         vjust = c(0, # Nordjlland
                   1, # Vestjylland
                   0, # Oestjylland
                   0, # Sjaelland
                   1, # Bornholm
                   1) # Fune_Lolland
  )

# Set map extent
main_panel_xlim <- st_bbox(biowide_regions)[c(1,3)] * c(0.6, 1.15)
main_panel_ylim <- st_bbox(biowide_regions)[c(2,4)] * c(0.98, 1.0125)
main_panel_width <- main_panel_xlim[2] - main_panel_xlim[1]
main_panel_height <- main_panel_ylim[2] - main_panel_ylim[1]

# Plot map
(map_plot <- ggplot() + 
    geom_sf(aes(colour = region,
                fill = region), 
            linewidth = 0.5,
            data = biowide_regions) +
    geom_richtext(aes(x = x, 
                      y = y, 
                      label = paste0("<span style = font-size:20pt>", region, "</span>"),
                      colour = region,
                      hjust = hjust,
                      vjust = vjust), 
                  data = biowide_regions,
                   fill = NA,
                  label.color = NA
              ) +
    scale_colour_manual(values = c("#C575D9", # Nordjlland 
                                   "#B88AC5", # Vestjylland
                                   "#7D3E8C", # Oestjylland 
                                   "#62B5B4", # Sjaelland 
                                   "#0F403F", # Bornholm 
                                   "#3D8A88"  # Fune_Lolland 
                                  )) +  
    scale_fill_manual(values = c("#F5FAFA", # Nordjlland
                                 "#F5FAFA", # Vestjylland
                                 "#F9F5FA", # Oestjylland
                                 "#F9F5FA", # Sjaelland
                                 "#F5FAFA", # Bornholm
                                 "#F9F5FA"  # Fune_Lolland
                                 )) + 
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


# Make overview map to show extent in the global context. 

# Get rnatural earth shapes
coast_lines <- ne_countries(scale = "medium") %>%
  st_union()
denmark <- ne_countries(scale = "medium", country = "Denmark")

# Determine map boundaries in Europe Lambert EPSG:3035
europe_lims <- data.frame(
  x = c(-10,40),
  y = c(30,65)) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3035) %>%
  st_bbox()

(overview_map <- ggplot() +
  geom_sf(data = coast_lines, 
          color = "grey40",
          fill = "grey70") +
  geom_sf(data = denmark, fill = "#7D3E8C") +
  geom_sf(data = denmark %>%
            st_transform(3035) %>%
            st_buffer(100 * 10^3) %>%
            st_bbox() %>%
            st_as_sfc(),  
          col = "black",
          linewidth = 4,
          lineend = "round",
          fill = "NA") +
  coord_sf(xlim = europe_lims[c(1,3)], 
           ylim = europe_lims[c(2,4)],
           crs = 3035) +
  theme_void())

save_plot("docs/figures/overview_map_figure_1.png",
          overview_map,
          base_asp = (europe_lims[3]-europe_lims[1]) / (europe_lims[4]-europe_lims[2]),
          bg = "white")

