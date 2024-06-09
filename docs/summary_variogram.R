# Summary variogram figure
# Jakob Assmann jakob.assmann@uzh.ch 6 June 2024

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)

# Load variogram files
load("data/variograms/variogram_list.Rda")
vario_df <- vario_list %>%
  map(function(x) x[[1]]) %>%
  bind_rows() %>%
  # filter(id %in%
  # c("amplitude_mean", 
  # "amplitude_sd", 	
  # "canopy_height",  
  # "Clay_utm32_10m", 
  # "dtm_10m", 	 		
  # "normalized_z_sd", 	
  # "ns_groundwater_summer_sd_110m",  
  # "openness_difference", 
  # "Sand_utm32_10m", 
  # "slope", 
  # "Soc_utm32_10m" ,
  # "solar_radiation",
  # "vegetation_density")) %>%
  split(., .$id) %>%
  map(function(x) mutate(x, gamma_std = scale(gamma, center = FALSE))) %>%
  bind_rows()

(ggplot(vario_df) +
  geom_line(aes(x = dist, y= gamma_std, colour = id)) +
    geom_vline(xintercept = 100, linetype = "dashed") +
    annotate("text", x = 100, y = Inf, hjust = 0, vjust = 1.5, 
             label = "  100 m") +
  labs(x = "Lag distance (m)", y = "Standardised semivariance",
       colour = "Predictor Variable") +
  theme_cowplot()) %>%
  save_plot("docs/variograms/st_variogram_1km.png", 
            .,
            base_asp = 2,
            bg = "white")
      