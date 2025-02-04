#!/usr/bin/env Rscript

# plot_map_and_pca.R
# 1) Reads geographic data + background map => "geo_plot" (colored by region).
# 2) Reads PCA data => "pca_plot_12" (PC1 vs PC2) and "pca_plot_34" (PC3 vs PC4).
# 3) Uses the same region color palette for all plots.
# 4) Arranges the two PCA plots (top row) above the geographic map (bottom row).
# 5) Saves final figure.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(sf)
  library(ggrepel)
  library(ggspatial)
  library(patchwork)   # For arranging plots
  library(stringr)     # For text normalization
})

# ─────────────────────────────────────────────────────────────────────────────
# 1) File Paths (edit as needed)
# ─────────────────────────────────────────────────────────────────────────────

BasePath    <- "~/Documents/genome_simulations_tsrelatedness/misc/geo_map/"
bg_map_file <- paste0(BasePath, "background_map_Jan2025_tsrel.rds")
pca_csv     <- paste0(BasePath, "balsac_branch_pca.csv")
output_file <- paste0(BasePath, "map_and_pca_grid2.jpg")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Define region color palette (shared)
# ─────────────────────────────────────────────────────────────────────────────

region_colors <- c(
  "Chaleur Bay"   = "#E41A1C",
  "Batiscan"      = "#377EB8",
  "Chaudière"     = "#4DAF4A",
  "L’Assomption"  = "#984EA3",   # Correct curly apostrophe
  "Mistassini"    = "#FF7F00"
)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Geographic Map Plot
# ─────────────────────────────────────────────────────────────────────────────

# Read background map
bg_map <- readRDS(bg_map_file)

# Hardcoded region coordinates
coords_means <- data.frame(
  X = c(-72.40168, -64.26971, -70.63634, -73.58594, -72.50327),
  Y = c(46.55456, 48.52418, 46.07026, 45.92777, 48.76101),
  region = c("Batiscan", "Chaleur Bay", "Chaudière", "L’Assomption", "Mistassini")
)

# Manually defined label positions for geographic regions
geo_labels <- data.frame(
  region = c("Chaleur Bay", "Batiscan", "Chaudière", "L’Assomption", "Mistassini"),
  X = c(-64.5, -72.1, -70.9, -73.2, -73.1),  # Adjust these manually
  Y = c(48.5, 47.1, 46.2, 45.2, 49.0)        # Adjust these manually
)

# Define a bounding box for easy manual adjustment
bbox <- st_bbox(c(
  xmin = -74.5,  # Western boundary
  xmax = -64.5,  # Eastern boundary
  ymin = 45.5,   # Southern boundary
  ymax = 49.3    # Northern boundary
), crs = 4326)

# Build the map plot
set.seed(123)  # Ensure consistent label positioning

geo_plot <- bg_map +
  geom_point(
    data = coords_means,
    aes(x = X, y = Y, color = region),
    size = 3
  ) +

  # Semi-transparent region label background
  geom_label_repel(
    data = geo_labels,
    aes(x = X, y = Y, label = region, color = region),
    fill = "white",
    force = 50,
    alpha = 0.6,
    color = NA,
    box.padding = 0.3,
    point.padding = 0.3,
    label.size = NA,
    segment.color = NA,
    seed = 123,
    size = 6
  ) +
  # Opaque text on top
  geom_label_repel(
    data = geo_labels,
    aes(x = X, y = Y, label = region, color = region),
    fill = NA,
    force = 50,
    box.padding = 0.3,
    point.padding = 0.3,
    label.size = NA,
    segment.color = NA,
    seed = 123,
    size = 6
  ) +

  coord_sf(
    xlim = c(bbox$xmin, bbox$xmax),
    ylim = c(bbox$ymin, bbox$ymax), 
    crs = 4326
  ) +

  # Scale bar & north arrow (bottom-right)
  annotation_scale(location = "br", width_hint = 0.3) +  # Move scale to bottom-left
  annotation_north_arrow(
    location = "br",  # Move north arrow to bottom-left
    which_north = "true",
    pad_x = unit(1.8, "in"), pad_y = unit(0.35, "in"),
    style = north_arrow_fancy_orienteering
  )+

  scale_color_manual(values = region_colors, guide = "none") +
  theme_bw() +
  theme(
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    plot.margin      = margin(0, 0, 0, 0),
    axis.line        = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "aliceblue")
  )

# Load required packages for the globe projection
suppressPackageStartupMessages({
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggspatial)
})

# Get a world map
world_map <- ne_countries(scale = "medium", returnclass = "sf")

# Transform world map to Lambert Azimuthal Equal-Area (LAEA) projection before cropping
world_map_laea <- st_transform(world_map, crs = "+proj=laea +lat_0=50 +lon_0=-70")

# Define bounding box for North America in LAEA projection
north_america_bbox <- st_bbox(c(
  xmin = -5000000, xmax = 2000000,  # Approximate extent in LAEA
  ymin = -2000000, ymax = 4000000
), crs = st_crs(world_map_laea))

# Crop the transformed world map
world_map_cropped <- st_crop(world_map_laea, north_america_bbox)

# Define the zoomed-in area as a proper polygon and transform it to LAEA
zoom_bbox <- st_as_sf(st_sfc(st_polygon(list(rbind(
  c(-74.5, 45.5),  # Lower-left
  c(-74.5, 49.3),  # Upper-left
  c(-64.5, 49.3),  # Upper-right
  c(-64.5, 45.5),  # Lower-right
  c(-74.5, 45.5)   # Close the polygon
))), crs = 4326)) %>%
  st_transform(crs = st_crs(world_map_laea))  # Transform to LAEA

# Create a properly zoomed-in inset map with consistent projection
inset_map <- ggplot() +
  geom_sf(
    data = world_map_cropped,
    fill = "lightgray", color = NA, size = 0.2  # Remove country borders
  ) +
  geom_sf(data = zoom_bbox, fill = "red", alpha = 0.3) +  # Highlight zoomed-in region
  coord_sf(crs = st_crs(world_map_laea), expand = FALSE) +  # Use same LAEA projection
  theme_void() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    plot.background = element_rect(fill = "aliceblue", color = "black")
  )

# Combine main map with inset
geo_plot_with_inset <- geo_plot +
  annotation_custom(
    grob = ggplotGrob(inset_map),
    xmin = bbox$xmax - 2,  # Move further right
    xmax = bbox$xmax,
    ymin = bbox$ymin,  # Push further down
    ymax = bbox$ymin + 1.5   # Maintain square shape
  )


# ─────────────────────────────────────────────────────────────────────────────
# 4) PCA Plots
# ─────────────────────────────────────────────────────────────────────────────

# Read PCA data and standardize region names
pca_df <- fread(pca_csv, header = TRUE) %>%
  mutate(region = str_replace_all(region, "'", "’"))  # Replace ASCII apostrophe

# Rename first four PCs
colnames(pca_df)[1:4] <- c("PC1", "PC2", "PC3", "PC4")

# Manually defined label positions for PCA plots
pca_labels <- data.frame(
  region = c("Chaleur Bay", "Batiscan", "Chaudière", "L’Assomption", "Mistassini"),
  PC1 = c(0.02, 0.01, 0.07, NA, NA),  # Adjust values manually for PC1 vs PC2
  PC2 = c(-0.07, 0.05, 0.03, NA, NA),
  PC3 = c(-0.06, 0.03, -0.05, 0.04, 0.09),  # Adjust values manually for PC3 vs PC4
  PC4 = c(0.05, 0.05, -0.035, -0.055, -0.01)
)

# PC1 vs. PC2 plot with manually placed labels
pca_plot_12 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(
    data = pca_labels,
    aes(x = PC1, y = PC2, label = region, color = region),
    size = 5
  ) +
  scale_color_manual(values = region_colors) +
  labs(x = "PC1", y = "PC2") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14)
  )

# PC3 vs. PC4 plot with manually placed labels
pca_plot_34 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = region)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(
    data = pca_labels,
    aes(x = PC3, y = PC4, label = region, color = region),
    size = 5
  ) +
  scale_color_manual(values = region_colors) +
  labs(x = "PC3", y = "PC4") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Combine Plots and Save
# ─────────────────────────────────────────────────────────────────────────────

# Arrange in a grid:
# Top row: PC1 vs PC2 | PC3 vs PC4
# Bottom row: Geo Map
combined_plot <- (pca_plot_12 | pca_plot_34) / geo_plot_with_inset +
  plot_layout(heights = c(1, 1)) + 
  plot_annotation(tag_levels = "A")

# Save final combined figure
ggsave(filename = output_file, plot = combined_plot, width = 10, height = 9, dpi = 300)
cat("Combined plot saved to:", output_file, "\n")