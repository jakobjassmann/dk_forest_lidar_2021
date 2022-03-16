# Figure 2 for the DK Forest LiDAR manuscript
# Jakob J. Assmann j.assmann@bio.au.dk 16 March 2022

# Dependencies
library(tidyverse)
library(ggplot)
library(patchwork)

# define patchwork plot layout
plot_layout <- "
AAAB
AAAC
AAAD
EFGH
"

# Generate map of high quality forests

  
# Set areas of interest
husby_klit <- st_bbox(445792.4, 6238762.3, 447477.4, 6240058.2) 