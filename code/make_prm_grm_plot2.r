#!/usr/bin/env Rscript

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(tools)
library(tibble)
library(scales)

set.seed(345)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: script.R grm_file.csv prm_file.csv metadata.csv")
}

# File paths
grm_file <- args[1]
prm_file <- args[2]
meta_file <- args[3]
grm_filename <- file_path_sans_ext(basename(grm_file))
prm_filename <- file_path_sans_ext(basename(prm_file))
output_dir <- dirname(grm_file)
plot_file <- file.path(output_dir, "grm_prm_heatmaps.jpg")

# ------------------------------------------------------------------------------
# 1) Read a symmetric NxN matrix from CSV
# ------------------------------------------------------------------------------
read_symmetric_matrix <- function(path) {
  df <- read_csv(path, col_names = TRUE)
  mat <- as.matrix(df)
  if (nrow(mat) == ncol(mat)) {
    rownames(mat) <- colnames(mat)
  } else {
    cat("WARNING: Matrix is not square! Dimensions:", nrow(mat), "x", ncol(mat), "\n")
  }
  mat
}

# ------------------------------------------------------------------------------
# 2) Process input files: return long data & clustered ranking
# ------------------------------------------------------------------------------
process_relatedness <- function(data_file, meta_file) {
  mat <- read_symmetric_matrix(data_file)
  metadata <- read_csv(meta_file, col_names = TRUE) %>%
    rename(depth = average_depth, proband = noid) %>%
    mutate(proband = as.character(proband)) %>%
    distinct(proband, proband_region, depth)
  matrix_ids <- rownames(mat)
  metadata <- metadata %>% filter(proband %in% matrix_ids)
  missing_in_meta <- setdiff(matrix_ids, metadata$proband)
  if (length(missing_in_meta) > 0) {
    cat("WARNING: These IDs are in the matrix but not in metadata:\n")
    print(missing_in_meta)
  }
  cat("Number of NA (missing) entries in this matrix:", sum(is.na(mat)), "\n")
  
  dist_mat <- as.dist(1 - mat)
  hc <- hclust(dist_mat, method = "average")
  clustered_proband <- rownames(mat)[hc$order]
  
    # (f) Region-based ordering with custom region order, then cluster ordering
  region_map <- metadata %>%
    select(proband, proband_region) %>%
    distinct() %>%
    mutate(proband_region = factor(proband_region,
                                   levels = c("L'Assomption", "Batiscan", "Chaudière", "Mistassini", "Chaleur Bay")))
  
  region_priority <- region_map %>%
    mutate(cluster_order = match(proband, clustered_proband)) %>%
    arrange(proband_region, cluster_order)
  
  proband_ranking <- region_priority$proband
  
  relatedness_long <- mat %>%
    as.data.frame() %>%
    rownames_to_column(var = "proband1") %>%
    pivot_longer(cols = -proband1, names_to = "proband2", values_to = "Relatedness") %>%
    mutate(Relatedness = ifelse(proband1 == proband2, NA, Relatedness))
  
  relatedness_long <- relatedness_long %>%
    left_join(metadata, by = c("proband1" = "proband")) %>%
    rename(proband_region1 = proband_region, depth1 = depth) %>%
    left_join(metadata, by = c("proband2" = "proband")) %>%
    rename(proband_region2 = proband_region, depth2 = depth)
  
  list(
    data = relatedness_long,
    ranking = proband_ranking
  )
}

# ------------------------------------------------------------------------------
# 3) Apply a common ranking to a long-format data frame
# ------------------------------------------------------------------------------
apply_common_ranking <- function(relatedness_data, common_proband_ranking) {
  relatedness_data %>%
    mutate(
      proband1 = factor(proband1, levels = common_proband_ranking),
      proband2 = factor(proband2, levels = common_proband_ranking)
    )
}

# ------------------------------------------------------------------------------
# 5) Region bar plot (manual color scheme)
# ------------------------------------------------------------------------------
generate_region_bar <- function(data) {
  region_colors <- c(
    "L'Assomption" = "#984EA3",
    "Batiscan"     = "#377EB8",
    "Chaudière"    = "#4DAF4A",
    "Mistassini"   = "#FF7F00",
    "Chaleur Bay"  = "#E41A1C"
  )
  
  ggplot(
    distinct(data, proband1, proband_region1, depth1),
    aes(x = proband1, y = depth1, fill = factor(proband_region1, 
                                                levels = c("L'Assomption", "Batiscan", "Chaudière", "Mistassini", "Chaleur Bay")))
  ) +
    geom_col(color = NA) +
    scale_fill_manual(
      values = region_colors,
      name = "Region",
      na.value = "blue",
      drop = FALSE,
      limits = c("L'Assomption", "Batiscan", "Chaudière", "Mistassini", "Chaleur Bay")
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 1))
}


# ------------------------------------------------------------------------------
# 6) Main script execution
# ------------------------------------------------------------------------------
res_grm <- process_relatedness(grm_file, meta_file)
res_prm <- process_relatedness(prm_file, meta_file)
common_proband_ranking <- res_prm$ranking
relatedness_grm <- res_grm$data
relatedness_prm <- res_prm$data

relatedness_grm <- apply_common_ranking(relatedness_grm, common_proband_ranking)
relatedness_prm <- apply_common_ranking(relatedness_prm, common_proband_ranking)

# --- Add indices based on common ranking ---
relatedness_grm <- relatedness_grm %>%
  mutate(row_index = match(as.character(proband1), common_proband_ranking),
         col_index = match(as.character(proband2), common_proband_ranking))
relatedness_prm <- relatedness_prm %>%
  mutate(row_index = match(as.character(proband1), common_proband_ranking),
         col_index = match(as.character(proband2), common_proband_ranking))

# --- Filter data for triangles ---
# For eGRM, use upper triangle (row_index > col_index)
upper_egrm <- relatedness_grm %>% filter(row_index > col_index)
# For PRM, use lower triangle (row_index < col_index)
lower_prm  <- relatedness_prm %>% filter(row_index < col_index)

# --- Create a common base with transparent background ---
base_transparent <- list(
  scale_x_discrete(limits = common_proband_ranking),
  scale_y_discrete(limits = rev(common_proband_ranking)),
  theme_void(),
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )
)

# --- Create individual ggplot calls for each triangle ---
# eGRM plot (bidirectional color scale) for upper triangle
plot_egrm <- ggplot(upper_egrm, aes(x = proband1, y = proband2, fill = Relatedness)) +
  geom_tile() +
  scale_fill_gradient2(
    name = "Branch GRM",
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
    trans = scales::pseudo_log_trans(sigma = 10),
    limits = c(-5000, 5000),
    oob = squish,
    na.value = "black",
    breaks = c(-1000, -100, 0, 100, 1000),
  ) +
  base_transparent
plot_egrm <- plot_egrm + theme(plot.title = element_blank())

# PRM plot (log scale) for lower triangle
plot_prm <- ggplot(lower_prm, aes(x = proband1, y = proband2, fill = Relatedness + 1e-4)) +
  geom_tile() +
  scale_fill_gradient(
    name = "Pedigree GRM",
    low = "white", high = "black", na.value = "black", trans = "log",
    breaks = 2^seq(-12, 0, by = 2),
    labels = format(2^seq(-12, 0, by = 2), digits = 3, scientific = FALSE)
  ) +
  base_transparent +
  guides(fill = guide_colorbar())  # Force legend to show
plot_prm <- plot_prm + theme(plot.title = element_blank())

# --- Extract legends from each plot (both set to "right") ---
legend_egrm <- cowplot::get_legend(
  plot_egrm + theme(legend.position = "right", legend.box.margin = margin(0,0,0,0))
)
legend_prm <- cowplot::get_legend(
  plot_prm + theme(legend.position = "right", legend.box.margin = margin(0,0,0,0))
)

# --- Combine the two legends vertically (eGRM on top, PRM below) ---
combined_legends <- plot_grid(legend_egrm, legend_prm, ncol = 1, rel_heights = c(1, 1))

# --- Remove legends from individual plots for overlay ---
plot_egrm_noleg <- plot_egrm + theme(legend.position = "none")
plot_prm_noleg <- plot_prm + theme(legend.position = "none")

# --- Overlay the two plots ---
overlay_plot <- ggdraw() +
  draw_plot(plot_prm_noleg, 0, 0, 1, 1) +
  draw_plot(plot_egrm_noleg, 0, 0, 1, 1)

# --- Combine overlay with region bar using old alignment ---
overlay_panel <- plot_grid(
  overlay_plot,
  generate_region_bar(relatedness_grm) + theme(plot.margin = margin(0, 0, 0, 5)),
  ncol = 1,
  rel_heights = c(1, 0.1),
  align = "v",
  axis = "lr"
)

# --- Combine overlay panel with the combined legends on the right ---
final_overlay <- plot_grid(
  overlay_panel,
  combined_legends,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

# --- Save the final composite ---
ggsave(plot_file, plot = final_overlay, width = 10, height = 10, dpi = 300)
cat("Overlayed heatmaps with legends, title, and aligned region bar saved to:", plot_file, "\n")


#Rscript ../../code/make_prm_grm_plot2.r  /Users/luke/Documents/genome_simulations_tsrelatedness/misc/geo_map/chr3_grm_noid.csv  /Users/luke/Documents/genome_simulations_tsrelatedness/misc/geo_map/chr3_prm_noid.csv /Users/luke/Documents/genome_simulations_tsrelatedness/misc/geo_map/balsac_proband_meta_noid.csv