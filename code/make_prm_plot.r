#!/usr/bin/env Rscript

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(tools)

set.seed(345)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the required input files are provided
if (length(args) < 3) {
  stop("Usage: script.R grm_file.csv prm_file.csv metadata.csv")
}

# Path to the input files
grm_file <- args[1]
prm_file <- args[2]
meta_file <- args[3]

# Extract filenames without path and extension
grm_filename <- file_path_sans_ext(basename(grm_file))
prm_filename <- file_path_sans_ext(basename(prm_file))

# Determine output file base name and location
output_dir <- dirname(grm_file)
plot_file <- file.path(output_dir, "grm_prm_heatmaps.jpg")

# Function to process input files
process_relatedness <- function(data_file, meta_file) {
  relatedness_matrix <- read_csv(data_file, col_names = TRUE)
  metadata <- read_csv(meta_file, col_names = TRUE)
  
  metadata <- metadata %>%
    rename(depth = average_depth, proband = noid) %>%
    distinct(proband, proband_region, depth)
  
  relatedness_matrix <- relatedness_matrix %>%
    mutate(proband1 = colnames(relatedness_matrix))
  
  mat <- as.matrix(relatedness_matrix[, -ncol(relatedness_matrix)])
  rownames(mat) <- relatedness_matrix$proband1
  dist_mat <- as.dist(1 - mat)
  hc <- hclust(dist_mat, method = "average")
  
  clustered_proband <- rownames(mat)[hc$order]
  region_map <- metadata %>% distinct(proband, proband_region) %>% filter(!is.na(proband_region))
  region_priority <- region_map %>% mutate(cluster_order = match(proband, clustered_proband)) %>%
    arrange(proband_region, cluster_order)
  
  proband_ranking <- region_priority$proband
  
  relatedness_long <- relatedness_matrix %>%
    pivot_longer(cols = -proband1, names_to = "proband2", values_to = "Relatedness") %>%
    mutate(Relatedness = ifelse(proband1 == proband2, NA, Relatedness),
           proband1 = as.numeric(proband1),
           proband2 = as.numeric(proband2))
  
  relatedness_long <- relatedness_long %>%
    left_join(metadata, by = c("proband1" = "proband")) %>%
    rename(proband_region1 = proband_region, depth1 = depth) %>%
    left_join(metadata, by = c("proband2" = "proband")) %>%
    rename(proband_region2 = proband_region, depth2 = depth)
  
  relatedness_long$proband1 <- factor(relatedness_long$proband1, levels = proband_ranking)
  relatedness_long$proband2 <- factor(relatedness_long$proband2, levels = proband_ranking)
  
  return(relatedness_long)
}

# Process both GRM and PRM
relatedness_grm <- process_relatedness(grm_file, meta_file)
relatedness_prm <- process_relatedness(prm_file, meta_file)

# Generate heatmaps
generate_heatmap <- function(data, title) {
  ggplot(data, aes(x = proband1, y = proband2, fill = Relatedness + 1e-4)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", na.value = 'blue', trans = "log") +
    theme_minimal() +
    labs(title = title, x = "Individuals", y = "Individuals") +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5))
}

heatmap_grm <- generate_heatmap(relatedness_grm, "Genetic Relatedness Matrix (GRM)")
heatmap_prm <- generate_heatmap(relatedness_prm, "Pedigree Relatedness Matrix (PRM)")

# Generate region color bars
region_colors_grm <- ggplot(
  distinct(relatedness_grm, proband1, proband_region1, depth1),
  aes(x = proband1, y = depth1, fill = proband_region1)
) +
  geom_col(color = NA) +
  scale_fill_viridis_d(option = "turbo", name = "Region", na.value = "blue") +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1))

region_colors_prm <- ggplot(
  distinct(relatedness_prm, proband1, proband_region1, depth1),
  aes(x = proband1, y = depth1, fill = proband_region1)
) +
  geom_col(color = NA) +
  scale_fill_viridis_d(option = "turbo", name = "Region", na.value = "blue") +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1))

# Combine heatmaps and their corresponding color bars
heatmap_grm_combined <- plot_grid(heatmap_grm, region_colors_grm, ncol = 1, rel_heights = c(1, 0.1), align = "v", axis = "lr")
heatmap_prm_combined <- plot_grid(heatmap_prm, region_colors_prm, ncol = 1, rel_heights = c(1, 0.1), align = "v", axis = "lr")

# Combine both heatmap-color bar pairs into final figure
combined_plot <- plot_grid(heatmap_grm_combined, heatmap_prm_combined, ncol = 2, align = "hv")

ggsave(plot_file, plot = combined_plot, width = 12, height = 7, dpi = 300)

cat("Heatmaps with aligned color bars saved to:", plot_file, "\n")
