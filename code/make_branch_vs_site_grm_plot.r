#!/usr/bin/env Rscript

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
library(tools)
library(tibble)

set.seed(345)

# ----------------------------------------------------------------------------
# Usage and arguments
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: script.R branch_grm.csv prm_grm.csv site_grm.csv metadata.csv")
}
branch_file <- args[1]  # Branch GRM (e.g., grm_branch_recap_noid.csv)
prm_file    <- args[2]  # Pedigree GRM (e.g., chr3_prm_noid.csv)
site_file   <- args[3]  # Site GRM (e.g., grm_site_recap_noid.csv)
meta_file   <- args[4]  # Metadata CSV (e.g., balsac_proband_meta_noid.csv)

# Output files
output_dir    <- dirname(branch_file)
overlay1_file <- file.path(output_dir, "grm_prm_heatmaps.jpg")
overlay2_file <- file.path(output_dir, "site_branch_grm_heatmaps.jpg")

# ----------------------------------------------------------------------------
# Function to read symmetric matrix
# ----------------------------------------------------------------------------
read_symmetric_matrix <- function(path) {
  df <- read_csv(path, col_names = TRUE)
  mat <- as.matrix(df)
  if (nrow(mat) != ncol(mat)) {
    stop("Matrix not square: ", path)
  }
  rownames(mat) <- colnames(mat)
  mat
}

# ----------------------------------------------------------------------------
# Process relatedness: hierarchical + region-based ranking + long data
# ----------------------------------------------------------------------------
process_relatedness <- function(mat_file, metadata_df) {
  mat <- read_symmetric_matrix(mat_file)
  meta <- metadata_df %>%
    rename(depth = average_depth, proband = noid) %>%
    mutate(proband = as.character(proband)) %>%
    distinct(proband, proband_region, depth)
  ids <- rownames(mat)
  meta <- filter(meta, proband %in% ids)

  # Hierarchical clustering
  dist_mat <- as.dist(1 - mat)
  hc <- hclust(dist_mat, method = "average")
  clustered <- rownames(mat)[hc$order]

  # Region-based ordering
  region_map <- meta %>%
    select(proband, proband_region) %>% distinct() %>%
    mutate(proband_region = factor(proband_region,
      levels = c("L'Assomption","Batiscan","Chaudière","Mistassini","Chaleur Bay")))
  region_priority <- region_map %>%
    mutate(cluster_order = match(proband, clustered)) %>%
    arrange(proband_region, cluster_order)
  ranking <- region_priority$proband

  # Long-format data
  long_df <- as.data.frame(mat) %>%
    rownames_to_column("p1") %>%
    pivot_longer(-p1, names_to = "p2", values_to = "Relatedness") %>%
    mutate(Relatedness = ifelse(p1 == p2, NA, Relatedness))

  list(data = long_df, ranking = ranking)
}

# ----------------------------------------------------------------------------
# Apply ranking to long df
# ----------------------------------------------------------------------------
apply_ranking <- function(df_long, ranking) {
  df_long %>% mutate(
    p1 = factor(p1, levels = ranking),
    p2 = factor(p2, levels = ranking)
  )
}

# ----------------------------------------------------------------------------
# Generate region-depth bar under heatmaps
# ----------------------------------------------------------------------------
generate_region_bar <- function(ranking, metadata_df) {
  meta <- metadata_df %>%
    rename(depth = average_depth, proband = noid) %>%
    mutate(proband = as.character(proband)) %>%
    distinct(proband, proband_region, depth)
  region_colors <- c(
    "L'Assomption" = "#984EA3",
    "Batiscan"     = "#377EB8",
    "Chaudière"    = "#4DAF4A",
    "Mistassini"   = "#FF7F00",
    "Chaleur Bay"  = "#E41A1C"
  )
  meta %>%
    filter(proband %in% ranking) %>%
    mutate(proband = factor(proband, levels = ranking)) %>%
    ggplot(aes(x = proband, y = depth, fill = proband_region)) +
    geom_col(color = NA) +
    scale_fill_manual(values = region_colors, limits = names(region_colors)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 1))
}

# ----------------------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------------------
# Read metadata
meta_df <- read_csv(meta_file, col_names = TRUE)

# Get ranking from pedigree GRM (prm_file)
res_prm <- process_relatedness(prm_file, meta_df)
ranking <- res_prm$ranking

# Process branch and site GRMs
res_branch <- process_relatedness(branch_file, meta_df)
res_site   <- process_relatedness(site_file,   meta_df)

df_branch <- apply_ranking(res_branch$data, ranking) %>% mutate(i = as.integer(p1), j = as.integer(p2))
df_prm    <- apply_ranking(res_prm$data,    ranking) %>% mutate(i = as.integer(p1), j = as.integer(p2))
df_site   <- apply_ranking(res_site$data,   ranking) %>% mutate(i = as.integer(p1), j = as.integer(p2))

# Define common theme
common_theme <- list(
  scale_x_discrete(limits = ranking),
  scale_y_discrete(limits = rev(ranking)),
  theme_void(),
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )
)

# Prepare triangles
upper_branch <- filter(df_branch, i > j)
lower_prm    <- filter(df_prm,    i < j)
lower_site   <- filter(df_site,   i < j)

# Plot 1: Branch (upper) + PRM (lower)
p_branch <- ggplot(upper_branch, aes(p1, p2, fill = Relatedness)) +
  geom_tile() +
  scale_fill_gradient2(
    name = "Branch GRM",
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
    trans = pseudo_log_trans(sigma = 10), limits = c(-5000,5000), oob = squish,
    na.value = "black"
  ) + common_theme + theme(plot.title = element_blank())

p_prm <- ggplot(lower_prm, aes(p1, p2, fill = Relatedness + 1e-4)) +
  geom_tile() +
  scale_fill_gradient(
    name = "Pedigree GRM",
    low = "white", high = "black", trans = "log",
    breaks = 2^seq(-12,0,2), labels = 2^seq(-12,0,2), na.value = "black"
  ) + common_theme + guides(fill = guide_colorbar()) + theme(plot.title = element_blank())

region_bar <- generate_region_bar(ranking, meta_df)

leg1 <- plot_grid(
  get_legend(p_branch + theme(legend.position = "right")),
  get_legend(p_prm   + theme(legend.position = "right")),
  ncol = 1
)

overlay1 <- ggdraw() +
  draw_plot(p_prm + theme(legend.position = "none"), 0,0,1,1) +
  draw_plot(p_branch + theme(legend.position = "none"), 0,0,1,1)
panel1 <- plot_grid(
  overlay1,
  region_bar,
  ncol = 1, rel_heights = c(1, 0.1)
)
final1 <- plot_grid(panel1, leg1, ncol = 2, rel_widths = c(1,0.2))

ggsave(overlay1_file, final1, width = 10, height = 10, dpi = 300)
cat("Saved Branch vs Pedigree overlay to:", overlay1_file, "\n")

# Plot 2: Branch (upper) + Site (lower)
p_site <- ggplot(lower_site, aes(p1, p2, fill = Relatedness + 1e-4)) +
  geom_tile() +
  scale_fill_gradient(
    name = "Site GRM",
    low = "#1ABC9C", high = "#E67E22", trans = "log",
    breaks = 2^seq(-12,0,2), labels = 2^seq(-12,0,2), na.value = "black"
  ) + common_theme + guides(fill = guide_colorbar()) + theme(plot.title = element_blank())

leg2 <- plot_grid(
  get_legend(p_branch + theme(legend.position = "right")),
  get_legend(p_site   + theme(legend.position = "right")),
  ncol = 1
)

overlay2 <- ggdraw() +
  draw_plot(p_site   + theme(legend.position = "none"), 0,0,1,1) +
  draw_plot(p_branch + theme(legend.position = "none"), 0,0,1,1)
panel2 <- plot_grid(
  overlay2,
  region_bar,
  ncol = 1, rel_heights = c(1, 0.1)
)
final2 <- plot_grid(panel2, leg2, ncol = 2, rel_widths = c(1,0.2))

ggsave(overlay2_file, final2, width = 10, height = 10, dpi = 300)
cat("Saved Branch vs Site overlay to:", overlay2_file, "\n")
