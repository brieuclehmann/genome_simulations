#!/usr/bin/env Rscript

# plot_deep_structure_towns.R
# 1) Reads "ts-relatedness-towns.csv" (already created in the first script).
# 2) Reads background map "background_map_April2022.rds".
# 3) Transforms coordinates and plots them with labels.
# 4) Saves the final figure "deep_structure_map.jpg".

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(sf)
  library(ggrepel)
  library(ggspatial)
})

# Hard-coded file paths (edit as needed)
towns_csv   <- "/Users/luke/Desktop/ts-relatedness-towns.csv"
bg_map_file <- "/Users/luke/Documents/Genizon/Data/RDS/background_map_Jan2025.rds"
output_file <- "~/Desktop/deep_structure_map.jpg"

# 1) Read the CSV and background map
towns_df <- fread(towns_csv)  # columns: region, name, t, lieu, Lat, Lon
bg_map   <- readRDS(bg_map_file)

# 2) Convert to sf and transform
crs_string <- "+proj=omerc +lat_0=46.8560266 +lonc=-71.6218555 +alpha=0 +k_0=.7 +datum=WGS84 +units=m +no_defs +gamma=35"

sf_towns <- towns_df %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4269) %>%
  st_transform(crs = crs_string)

# Convert to data frame for ggplot
coords_df <- as.data.frame(st_coordinates(sf_towns))
coords_df$region <- towns_df$region
coords_df$name   <- towns_df$name
coords_df$t      <- towns_df$t
coords_df$lieu   <- towns_df$lieu

# We'll build a label, e.g. "region name"
coords_df$LAB <- paste(coords_df$region, coords_df$name, coords_df$lieu)

# 3) Construct the plot
geo_plot <- bg_map +
  geom_point(data = coords_df, aes(x = X, y = Y), size = 0.3, color = "black") +
  geom_text_repel(
    data = coords_df,
    aes(x = X, y = Y, label = LAB),
    size = 1.5,
    max.overlaps = Inf,
    force = 5,
    nudge_y = 0.5,
    direction = "both",
    min.segment.length = 0,
    segment.size = 0.2,
    box.padding = 0.2,
    point.padding = 0.3
  ) +
  coord_sf(
    crs = crs_string,
    xlim = c(-249000, 382000),
    ylim = c(-118000, 160000)
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.title= element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.title    = element_blank(),
    legend.position = "bottom",
    legend.direction= "horizontal",
    legend.text     = element_text(size = rel(0.8)),
    plot.margin     = margin(0, 0, 0, 0),
    axis.line       = element_blank(),
    panel.grid.major= element_blank(),
    panel.background= element_rect(fill = "aliceblue")
  ) +
  ggtitle("Deep Structure Parishes")

# 4) Save final figure
ggsave(geo_plot, filename = output_file, width = 7, height = 4)
cat("Plot saved to:", output_file, "\n")
