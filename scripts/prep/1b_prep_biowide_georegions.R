library(sf)
library(tidyverse)
library(ggplot2)
setwd("data/biowide_georegions/")
biowide_geoms <- list.files(pattern = "shp$", recursive = T) %>%
  map(read_sf) %>%
  map(st_transform, crs = 25832) %>%
  set_names(c("DK", "Nordjylland", "Oestjylland", "Sjaelland_clip", "Vestjylland"))
denmark <- pluck(biowide_geoms, 1)
regions <- biowide_geoms[2:length(biowide_geoms)]
plot(st_geometry(denmark))
plot(regions[[]])
ggplot() +
  geom_sf(data = denmark) +
  geom_sf(data = regions[[1]], fill = "red") +
  geom_sf(data = regions[[2]], fill = "green") +
  geom_sf(data = regions[[3]], fill = "blue") +
  geom_sf(data = regions[[4]], fill = "yellow")

# Set a buffer to add some tolerance to the boundaries
buffer <- 100
# Clean up geometries
Nordjlland <- biowide_geoms$DK %>%
  st_buffer(buffer) %>%
  st_crop(st_bbox(biowide_geoms$Nordjylland) + c(-buffer, 0, buffer, buffer)) %>%
  st_union() %>%
  st_sf() %>%
  mutate(region = "Nordjlland")
plot(st_geometry(Nordjlland))
Vestjylland <- biowide_geoms$DK %>%
  st_buffer(buffer) %>%
  st_crop(st_bbox(biowide_geoms$Vestjylland) + c(-buffer, -buffer, 0, 0)) %>%
  st_union() %>%
  st_sf() %>%
  mutate(region = "Vestjylland")

plot(st_geometry(Vestjylland))
# Build Oestjylland using a manual cut out based on original Oestylland geoms
Jutland <- biowide_geoms$DK %>%
  filter(OBJECTID == 159) %>%
plot(st_geometry(Jutland))
Oestjylland_mask <- st_polygon(list(
  matrix(
    c(595531.160,6044625.989,
      549985.717,6127215.646,
      541352.345,6134615.679,
      544479.740,6146816.924,
      542365.445,6148622.885,
      542497.588,6153159.810,
      547144.633,6152190.758,
      547056.537,6156265.181,
      549236.904,6156375.300,
      575489.403,6172981.327,
      608260.979,6174214.666,
      682613.694,6291734.242,
      st_bbox(Nordjlland)["xmax"], st_bbox(Nordjlland)["ymin"],
      st_bbox(Vestjylland)["xmax"], st_bbox(Nordjlland)["ymin"],
      st_bbox(Vestjylland)["xmax"], st_bbox(Jutland)["ymin"] - buffer,
      595531.160,6044625.989),
    ncol = 2,
    byrow = T)))
Oestjylland <-  Oestjylland_mask %>%
  st_intersects(biowide_geoms$DK) %>%
    unlist() %>%
    biowide_geoms$DK[.,] %>%
    st_buffer(buffer) %>%
    st_crop(xmin = as.numeric(st_bbox(Oestjylland_mask)["xmin"]),
            ymin = as.numeric(st_bbox(Oestjylland_mask)["ymin"]),
            xmax = as.numeric(st_bbox(Oestjylland_mask)["xmax"]),
            ymax = as.numeric(st_bbox(Nordjlland)["ymin"])) %>%
  st_union() %>%
  st_sf() %>%
  mutate(region = "Oestjylland")

# Sjaelland
Sjaelland <- biowide_geoms$Sjaelland_clip %>%
  st_intersects(biowide_geoms$DK) %>%
  unlist() %>%
  biowide_geoms$DK[.,] %>%
  filter(!(OBJECTID %in% c(1521, 1547, 1544, 1549, 1555, 1510))) %>%
  st_buffer(buffer) %>%
  st_union() %>%
  st_sf() %>%
  mutate(region = "Sjaelland")

# Bornholm
Bornholm <- biowide_geoms$DK %>%
  st_crop(xmin = 830383.108,
          ymin = as.numeric(st_bbox(biowide_geoms$DK)["ymin"]),
          xmax = as.numeric(st_bbox(biowide_geoms$DK)["xmax"]),
          ymax = as.numeric(st_bbox(biowide_geoms$DK)["ymax"])) %>%
  st_buffer(buffer) %>%
  st_union() %>%
  st_sf() %>%
  mutate(region = "Bornholm")

# Compile biowide zones
biowide_zones <- bind_rows(
  Nordjlland,
  Vestjylland,
  Oestjylland,
  Sjaelland,
  Bornholm)
plot(biowide_zones)

# Fune & Lolland (by exclusion)
fune_lolland <- biowide_zones %>%
  st_intersects(biowide_geoms$DK) %>%
  unlist() %>%
  map_dbl(~.x*-1) %>%
  biowide_geoms$DK[.,] %>%
  st_buffer(buffer) %>%
  st_union() %>%
  st_sf() %>%
  mutate(region = "Fune_Lolland")
plot(st_geometry(fune_lolland))

# Combine into one sf object and write out as shapefile
biowide_zones <- bind_rows(
  Nordjlland,
  Vestjylland,
  Oestjylland,
  Sjaelland,
  Bornholm,
  fune_lolland)
write_sf(biowide_zones, "biowide_zones.shp")
