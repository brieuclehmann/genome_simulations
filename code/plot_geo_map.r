#!/usr/bin/env Rscript

# plot_map_and_pca.R
# 1) Reads "ts-relatedness-towns.csv" for geographic map data
# 2) Reads "balsac_branch_pca.csv" for PCA data
# 3) Uses a shared color palette for region
# 4) Creates side-by-side plots: (A) Map, (B) PCA scatter
# 5) Saves final figure as "map_pca_combined.jpg"

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(sf)
  library(ggrepel)
  library(ggspatial)
  library(patchwork)  # For combining plots side-by-side
  # Alternatively: library(cowplot)
})

# ─────────────────────────────────────────────────────────────────────────────
# 1) File paths (edit as needed)
# ─────────────────────────────────────────────────────────────────────────────
towns_csv     <- "/Users/luke/Desktop/ts-relatedness-towns.csv"
bg_map_file   <- "/Users/luke/Documents/Genizon/Data/RDS/background_map_Jan2025.rds"
pca_csv       <- "/Users/luke/Desktop/tsrelatedness/balsac_branch_pca.csv"
output_file   <- "~/Desktop/map_pca_combined.jpg"

# ─────────────────────────────────────────────────────────────────────────────
# 2) Region-based color palette
#    Make sure each region in your data appears here.
# ─────────────────────────────────────────────────────────────────────────────
region_cols <- c(
  "Chaleur Bay"  = "#F8766D",
  "Batiscan"     = "#A3A500",
  "Chaudière"    = "#00BF7D",
  "L’Assomption" = "#00B0F6",
  "Mistassini"   = "#E76BF3"
)
# If your data includes more/different regions, expand or adjust region_cols.

# ─────────────────────────────────────────────────────────────────────────────
# 3) Read the map data (towns + background map)
# ─────────────────────────────────────────────────────────────────────────────
message("Reading geographic data...")
towns_df <- fread(towns_csv)   # columns: region, name, t, lieu, Lat, Lon
bg_map   <- readRDS(bg_map_file)  # a ggplot object with background

# Convert to sf and transform coordinates
crs_string <- "+proj=omerc +lat_0=46.8560266 +lonc=-71.6218555 +alpha=0 +k_0=.7 +datum=WGS84 +units=m +no_defs +gamma=35"

sf_towns <- towns_df %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4269) %>%  # adjust 4269 if your lat/lon is WGS84(4326)
  st_transform(crs = crs_string)

# Create a data frame for plotting
coords_df <- as.data.frame(st_coordinates(sf_towns))
coords_df$region <- towns_df$region
coords_df$name   <- towns_df$name

# ─────────────────────────────────────────────────────────────────────────────
# 4) Create the map plot
# ─────────────────────────────────────────────────────────────────────────────
map_plot <- bg_map +
  geom_point(
    data = coords_df,
    aes(x = X, y = Y, color = region),
    size = 2
  ) +
  coord_sf(
    crs = crs_string,
    xlim = c(-249000, 382000),
    ylim = c(-118000, 160000)
  ) +
  scale_color_manual(values = region_cols) +
  theme_bw() +
  theme(
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    plot.title      = element_text(hjust = 0.5),
    plot.margin     = margin(5,5,5,5),
    legend.title    = element_blank(),
    legend.position = "bottom",
    legend.direction= "horizontal",
    panel.grid.major= element_blank(),
    panel.background= element_rect(fill = "aliceblue")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Read PCA data
# ─────────────────────────────────────────────────────────────────────────────
# The CSV has 11 columns: PC1..PC10 + region
# We'll rename them explicitly:
message("Reading PCA data...")
pca_df <- read.csv(pca_csv, header = FALSE)
colnames(pca_df) <- c(
  "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","region"
)

print(head(pca_df))
# If your file has a header row, you can skip this rename or set col_names in read_csv.

# ─────────────────────────────────────────────────────────────────────────────
# 6) Create the PCA plot
# ─────────────────────────────────────────────────────────────────────────────
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = region_cols) +
  labs(x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.margin     = margin(5,5,5,5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 7) Combine the plots side-by-side using cowplot
# ─────────────────────────────────────────────────────────────────────────────
message("Combining map plot and PCA plot using cowplot...")

# Align the legend between the two plots
legend <- get_legend(pca_plot + theme(legend.position = "bottom"))

# Remove the legends from the individual plots
map_plot_clean <- map_plot + theme(legend.position = "none")
pca_plot_clean <- pca_plot + theme(legend.position = "none")

# Arrange map (smaller) and PCA (larger) with aligned legend below
combined_plot <- plot_grid(
  map_plot_clean, pca_plot_clean, legend,
  ncol = 2, rel_widths = c(0.4, 1),
  align = "hv"
)

# ─────────────────────────────────────────────────────────────────────────────
# 8) Save the final figure
# ─────────────────────────────────────────────────────────────────────────────
ggsave(output_file, combined_plot, width = 14, height = 6, dpi = 300)
message("Saved combined plot to: ", output_file)
