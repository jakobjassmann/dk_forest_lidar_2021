# Script for DK Forest LiDAR Figure 2
# Jakob J. Assmann j.assmann@bio.au.dk 23 March 2022

# Dependencies 

library(terra)
library(sf)
library(caret)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(cowplot)
library(ggnewscale)

# Set colours
high_quality_col <- "#3CB0AE"
low_quality_col <- "#D67D49" 

# Load pixel sample
load("data/training_data/pixel_training.Rda")

# Get Denmark boundaries
denmark <- read_sf("data/stratification/biowide_georegions/DK/DK.shp")

# Load biowide regions 
biowide_regions <- read_sf("data/stratification/biowide_georegions/biowide_zones.shp")

# Load validation data
load("data/training_data/pixel_valid_biowide.Rda")

# Load training data
load("data/training_data/pixel_training_biowide.Rda")

# Load best model
load("data/models/final_ranger_model_pixel_biowide.Rda")
rf_biowide <- rf_fit
rm(rf_fit)

# Calculate performances based on  validation data
performance_overall <- predict(rf_biowide, 
                                       newdata = pixel_valid_biowide) %>%
  confusionMatrix(data = ., 
                  pixel_valid_biowide$forest_value)

performance_list <- pixel_valid_biowide %>% 
  ungroup() %>%
  split(., .$biowide_region) %>%
  lapply(function(test_data_region){
    test_preds_region <- predict(rf_biowide, newdata = test_data_region)
    #cat(unique(test_data_region$biowide_region), "\n")
    confusionMat <- confusionMatrix(data = test_preds_region, test_data_region$forest_value)
    #print(confusionMat)
    return(confusionMat)
  }) 


# Extract statistics from performance_list
row_order <- sapply(biowide_regions$region, function (x) which(x == names(performance_list)))
biowide_regions$accuracy <- sapply(performance_list, 
                                   function(x) round(x$overall["Accuracy"],2))[row_order]
biowide_regions$sens <- sapply(performance_list, 
                               function(x) round(x$byClass["Sensitivity"],2))[row_order]


biowide_regions$useracc <- sapply(performance_list, 
                                  function(x) round(x$byClass["Pos Pred Value"],2))[row_order]

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
         x = x_cen + map_span_x * c(-0.25,  # Nordjlland
                                    -0.15,  # Vestjylland
                                    0.30,  # Oestjylland
                                    0.25,  # Sjaelland
                                    -0.10,  # Bornholm
                                    -0.30), # Fune_Lolland
         
         y = y_cen + map_span_y * c( -0.10,  # Nordjlland
                                     0.10,  # Vestjylland
                                     0.30,  # Oestjylland
                                     0.05,  # Sjaelland
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

# Update region names
biowide_regions$region <- factor(
  c("Northern<br>Jutland", 
    "Western<br>Jutland", 
    "Eastern<br>Jutland", 
    "Zealand", 
    "Bornholm", 
    "Funen-Lolland-Falster-Møn"),
  levels = 
    c("Northern<br>Jutland", 
      "Western<br>Jutland", 
      "Eastern<br>Jutland", 
      "Zealand", "Bornholm", 
      "Funen-Lolland-Falster-Møn")
)

# Generate labels
biowide_regions <- mutate(
  biowide_regions,
  label = paste0("<span style='font-size:20pt'>",
                 #"**", 
                 region, 
                 #"**",
                 "</span>",
                 "<span style='font-size:14pt; color:black'>",
                 "<br>Accuracy: ", formatC(accuracy, digits = 2, format = "f"),
                 "<br>Sensitivity: ", formatC(sens, digits = 2, format = "f"),
                 "<br>User Accuracy: ", formatC(useracc, digits = 2, format = "f"),
                 "</span>"))

# Build forest type panels

# Set known coordinates (selected manually)
p25 <- c(576102,  6219809)
old_growth <- c(576699,  6196176)
plantation <- c(462561, 6225310)

# Expand to 100 m x 100 m
p25 <- st_bbox(c(xmin = p25[1] - 64, xmax = p25[1] + 64,
                 ymin = p25[2] - 50, ymax = p25[2] + 50),
               crs = st_crs(biowide_regions)) 
old_growth <- st_bbox(c(xmin = old_growth[1] - 64, xmax = old_growth[1] + 64,
                        ymin = old_growth[2] - 50, ymax = old_growth[2] + 50),
                      crs = st_crs(biowide_regions)) 
plantation <- st_bbox(c(xmin = plantation[1] - 64, xmax = plantation[1] + 64,
                        ymin = plantation[2] - 50, ymax = plantation[2] + 50),
                      crs = st_crs(biowide_regions))

# Generate orthos 
# p25_ortho <- get_ortho(p25, "p25", scale = FALSE)
# old_growth_ortho <- get_ortho(old_growth, "old_growth", scale = FALSE)
# plantation_ortho <- get_ortho(plantation, "plantation", scale = FALSE)

# Load orthos
p25_ortho <- rast("data/orthophotos/p25.tif")
old_growth_ortho <- rast("data/orthophotos/old_growth.tif")
plantation_ortho <- rast("data/orthophotos/plantation.tif")

plot_ortho <- function(ortho, ortho_name, frame_colour = "black",
                       scale_caption = "50 m"){
  width <- ext(ortho)[2] - ext(ortho)[1]
  height <- ext(ortho)[4] - ext(ortho)[3]
  temp_file <- tempfile()
  png(temp_file, width = 640 * 2, height = 500 * 2)
  plotRGB(ortho)
  text(ext(ortho)[1] + width * 0.04,
       ext(ortho)[3] + height * 0.88,
       ortho_name,
       adj = 0,
       col = "white",
       cex = 10)
  rect(ext(ortho)[1] + width * 0.9 - 50,
       ext(ortho)[3] + height * 0.1,
       ext(ortho)[1] + width * 0.9,
       ext(ortho)[3] + height * 0.1 + height * 0.025,
       col = "white",
       border = "white")
  text(ext(ortho)[1] + width * 0.9 - 25,
       ext(ortho)[3] + height * 0.21,
       scale_caption,
       col = "white",
       cex = 10)
  dev.off()
  gg_grob <- ggplot() +
    draw_image(temp_file) +
    geom_rect(aes(xmin = 0, xmax = 1,
             ymin = 0.5 - (0.5 * (500) / (640)), 
             ymax = 0.5 + (0.5 * (500) / (640))), 
             colour = frame_colour,
             size = 1.5,
             fill = NA) +
     # labs(subtitle = ortho_name) +
    coord_equal() +
    theme_map() +
    theme(plot.margin=unit(c(0,0,0,0),"mm"),
          #panel.border = element_rect(colour = "red", fill = NA)
          )
  rm(temp_file)
  return(gg_grob)
}

p25_grob <- plot_ortho(p25_ortho, "§25 Forest", high_quality_col)
old_growth_grob <- plot_ortho(old_growth_ortho, "Private Untouched",
                              high_quality_col)
plantation_grob <- plot_ortho(plantation_ortho, "Plantation",
                              low_quality_col)

# Set map extent
main_panel_xlim <- st_bbox(biowide_regions)[c(1,3)] * c(0.70, 1.05)
main_panel_ylim <- st_bbox(biowide_regions)[c(2,4)] * c(0.99, 1.0025)
main_panel_width <- main_panel_xlim[2] - main_panel_xlim[1]
main_panel_height <- main_panel_ylim[2] - main_panel_ylim[1]

main_panel_height/main_panel_width

main_panel <- 
  ggplot() + 
  geom_sf(#aes(colour = region,
          #    fill = region), 
          size = 0.5,
          fill = NA,
          data = biowide_regions) +
  geom_richtext(aes(x = x, 
                    y = y, 
                    label = label,
                    #colour = region,
                    hjust = hjust,
                    vjust = vjust), 
                data = biowide_regions,
                fill = NA,
                label.color = NA) +
  geom_sf(data = pixel_training_data, # %>% sample_n(20), 
          mapping = aes(colour = forest_value),
          size = 0.001) +
  scale_colour_manual(values = c(high_quality_col,
                                 low_quality_col)) +
  coord_sf(clip = "off",
           crs = st_crs(biowide_regions),
           xlim = main_panel_xlim,
           ylim = main_panel_ylim) +
  labs(title = "Overall Performance - Best Model: Random Forest",
       subtitle = paste0(
         "Accuracy: ", round(performance_overall$overall["Accuracy"], 2), 
         "\nSensitivity: ", round(performance_overall$byClass["Sensitivity"], 2),
         "\nUser Accuracy: ", round(performance_overall$byClass["Pos Pred Value"], 2))) +
  annotate("text", 
           x = main_panel_xlim[1] + 
             0.63 * main_panel_width,
           y = 11000 + main_panel_ylim[1] +
             1.15 * main_panel_height,
           label = "Training: High Value Pixels", 
           colour = "black",
           size = 14 * 0.35,
           hjust = 0,
           vjust = 0.5) +
  annotate("rect", 
           xmin = main_panel_xlim[1] + 
             0.59 * main_panel_width,
           xmax = 20000 + main_panel_xlim[1] + 
             0.59 * main_panel_width,
           ymin = main_panel_ylim[1] +
             1.15 * main_panel_height,
           ymax = 20000 + main_panel_ylim[1] +
             1.15 * main_panel_height,
           color = "black",
           fill = high_quality_col) +
  annotate("text", 
           x = main_panel_xlim[1] + 
             0.63 * main_panel_width,
           y = 11000 + main_panel_ylim[1] +
             1.09 * main_panel_height,
           label = "Training: Low Value Pixels", 
           colour = "black",
           size = 14 * 0.35,
           hjust = 0,
           vjust = 0.5) +
  annotate("rect", 
           xmin = main_panel_xlim[1] + 
             0.59 * main_panel_width,
           xmax = 20000 + main_panel_xlim[1] + 
             0.59 * main_panel_width,
           ymin = main_panel_ylim[1] +
             1.09 * main_panel_height,
           ymax = 20000 +main_panel_ylim[1] +
             1.09 * main_panel_height,
           color = "black",
           fill = low_quality_col) +
  annotate("rect",
           xmin =main_panel_xlim[1] + 
             0.01 * main_panel_width,
           xmax = 100000 + main_panel_xlim[1] + 
             0.01 * main_panel_width,
           ymin = main_panel_ylim[1] + 5000,
           ymax = main_panel_ylim[1] + 5000 + 5000,
           fill = "black") +
  annotate("text", 
           x = 50000 + main_panel_xlim[1] + 
             0.01 * main_panel_width,
           y = main_panel_ylim[1] + 5000 + 25000,
           label = "100 km", 
           colour = "black",
           size = 14 * 0.35,
           hjust = 0.5,
           vjust = 0.5) +
  theme_map() +
  theme(legend.position = "none",
        plot.margin = margin(0.25,0,0.25,0.25, unit = "in"),
        #panel.border = element_rect(colour = "red", fill = NA)
        )

plot_grid(main_panel,
          plot_grid(p25_grob,
                    old_growth_grob,
                    plantation_grob,
                    nrow = 3),
          ncol = 2,
          rel_widths = c(2.5,1)) %>%
  save_plot("docs/figures/figure_2.png", 
          .,
          base_height = 6,
          bg = "white")

