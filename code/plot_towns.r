#!/usr/bin/env Rscript

# plot_deep_structure_towns.R
# 1) Reads "ts-relatedness-towns.csv" (already created in the first script).
# 2) Reads background map "background_map_Jan2025.rds".
# 3) Transforms coordinates and plots them with semi-transparent region labels.
# 4) Saves the final figure "deep_structure_map.jpg".

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(sf)
  library(ggrepel)
  library(ggspatial)
})

# ─────────────────────────────────────────────────────────────────────────────
# 1) File Paths (Modify as Needed)
# ─────────────────────────────────────────────────────────────────────────────
towns_csv   <- "~/Documents/genome_simulations_tsrelatedness/misc/geo_map/ts-relatedness-towns.csv"
bg_map_file <- "~/Documents/genome_simulations_tsrelatedness/misc/geo_map/background_map_Jan2025_tsrel.rds"
output_file <- "~/Documents/genome_simulations_tsrelatedness/misc/geo_map/deep_structure_map.jpg"

# ─────────────────────────────────────────────────────────────────────────────
# 2) Read the CSV and Background Map
# ─────────────────────────────────────────────────────────────────────────────
towns_df <- fread(towns_csv)  # columns: region, name, t, lieu, Lat, Lon
bg_map   <- readRDS(bg_map_file)

# Convert to sf and transform
crs_string <- "+proj=omerc +lat_0=46.8560266 +lonc=-71.6218555 +alpha=0 +k_0=.7 +datum=WGS84 +units=m +no_defs +gamma=35"

sf_towns <- towns_df %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4269) %>%
  st_transform(crs = crs_string)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Compute Region Centroids for Labels
# ─────────────────────────────────────────────────────────────────────────────
sf_region_centroids <- sf_towns %>%
  group_by(region) %>%
  summarise(geometry = st_centroid(st_union(geometry)), .groups = "drop")

# Convert centroid sf to data frame for ggplot labels
coords_region <- as.data.frame(st_coordinates(sf_region_centroids))
coords_region$region <- sf_region_centroids$region

# ─────────────────────────────────────────────────────────────────────────────
# 4) Construct the Plot
# ─────────────────────────────────────────────────────────────────────────────

geo_plot <- bg_map +
  geom_point(data = as.data.frame(st_coordinates(sf_towns)), aes(x = X, y = Y), 
             size = 0.3, color = "black") +
  
  # Semi-transparent label background (alpha = 0.4)
  geom_label_repel(
    data = coords_region,
    aes(x = X, y = Y, label = region),
    size = 5,
    fontface = "bold",
    fill = "white",
    alpha = 0.6, # Transparent background
    color = NA,  # Hide text on this layer
    box.padding = 0.3,
    point.padding = 0.3,
    label.size = NA, # Removes border
    segment.color = NA, # Removes segment line
    seed = 123 # Ensures consistent label positioning
  ) +

  # Opaque text labels with transparent background
  geom_label_repel(
    data = coords_region,
    aes(x = X, y = Y, label = region),
    size = 5,
    fontface = "bold",
    fill = NA,  # Transparent label background
    color = "black",  # Fully opaque text
    box.padding = 0.3,
    point.padding = 0.3,
    label.size = NA, # Removes border
    segment.color = NA, # Removes segment line
    seed = 123 # Ensures consistent label positioning
  ) +

  # Coordinate and theme settings
  coord_sf(
    crs = crs_string,
    xlim = c(-249000, 386000),
    ylim = c(-118000, 160000)
  ) +
  
  # Scale bar & north arrow (now bottom left)
  annotation_scale(
    location = "bl", # Bottom-left corner
    width_hint = 0.5
  ) +
  annotation_north_arrow(
    location = "bl", # Bottom-left corner
    which_north = "true",
    pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering
  ) +

  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.title= element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title    = element_blank(),
    legend.position = "bottom",
    legend.direction= "horizontal",
    legend.text     = element_text(size = rel(0.8)),
    plot.margin     = margin(0, 0, 0, 0),
    axis.line       = element_blank(),
    panel.grid.major= element_blank(),
    panel.background= element_rect(fill = "aliceblue")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Save Final Figure
# ─────────────────────────────────────────────────────────────────────────────
ggsave(geo_plot, filename = output_file, width = 8.5, height = 3.8, dpi = 300)
cat("Plot saved to:", output_file, "\n")
