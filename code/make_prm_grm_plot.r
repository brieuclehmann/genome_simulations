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


# ------------------------------------------------------------------------------
# 1) Function to read a symmetric NxN matrix from CSV
#    Row i => ID i, Col i => ID i, first row is column names
# ------------------------------------------------------------------------------
read_symmetric_matrix <- function(path) {
  # Read CSV, specifying col_names=TRUE so first line is column names (IDs)
  df <- read_csv(path, col_names = TRUE)
  
  # Convert to matrix
  mat <- as.matrix(df)
  
  # The matrix is NxN, with colnames the IDs; set rownames = colnames
  if (nrow(mat) == ncol(mat)) {
    rownames(mat) <- colnames(mat)
  } else {
    cat("WARNING: Matrix is not square! Dimensions:", 
        nrow(mat), "x", ncol(mat), "\n")
  }
  mat
}

# ------------------------------------------------------------------------------
# 2) Function to process input files: return both the long data and a clustered ranking
# ------------------------------------------------------------------------------
process_relatedness <- function(data_file, meta_file) {
  
  # (a) Read NxN matrix
  mat <- read_symmetric_matrix(data_file)
  
  # (b) Read and prepare metadata
  metadata <- read_csv(meta_file, col_names = TRUE) %>%
    rename(depth = average_depth, proband = noid) %>%
    mutate(proband = as.character(proband)) %>%
    distinct(proband, proband_region, depth)
  
  # (c) Subset metadata to only keep IDs that appear in the matrix
  matrix_ids <- rownames(mat)
  metadata <- metadata %>% filter(proband %in% matrix_ids)
  
  # (d) Sanity checks
  missing_in_meta <- setdiff(matrix_ids, metadata$proband)
  if (length(missing_in_meta) > 0) {
    cat("WARNING: These IDs are in the matrix but not in metadata:\n")
    print(missing_in_meta)
  }
  
  num_nas <- sum(is.na(mat))
  cat("Number of NA (missing) entries in this matrix:", num_nas, "\n")
  
  # (e) Build hierarchical clustering-based ranking
  dist_mat <- as.dist(1 - mat)
  hc <- hclust(dist_mat, method = "average")
  clustered_proband <- rownames(mat)[hc$order]

  # (f) Region-based ordering first, then cluster ordering
  region_map <- metadata %>%
    select(proband, proband_region) %>%
    distinct()
  
  region_priority <- region_map %>%
    mutate(cluster_order = match(proband, clustered_proband)) %>%
    arrange(proband_region, cluster_order)
  
  proband_ranking <- region_priority$proband
  
  # (g) Convert NxN matrix into long format
  relatedness_long <- mat %>%
    as.data.frame() %>%
    rownames_to_column(var = "proband1") %>%
    pivot_longer(
      cols = -proband1,
      names_to = "proband2",
      values_to = "Relatedness"
    ) %>%
    mutate(
      # Remove diagonal
      Relatedness = ifelse(proband1 == proband2, NA, Relatedness)
    )
  
  # (h) Merge with metadata (attach region/depth)
  relatedness_long <- relatedness_long %>%
    left_join(metadata, by = c("proband1" = "proband")) %>%
    rename(proband_region1 = proband_region, depth1 = depth) %>%
    left_join(metadata, by = c("proband2" = "proband")) %>%
    rename(proband_region2 = proband_region, depth2 = depth)
  
  # Return both
  list(
    data = relatedness_long,
    ranking = proband_ranking
  )
}

# ------------------------------------------------------------------------------
# 3) Helper function to apply a common ranking to a long-format data frame
# ------------------------------------------------------------------------------
apply_common_ranking <- function(relatedness_data, common_ranking) {
  relatedness_data %>%
    mutate(
      proband1 = factor(proband1, levels = common_ranking),
      proband2 = factor(proband2, levels = common_ranking)
    )
}

# ------------------------------------------------------------------------------
# 4) Heatmap generation function that allows linear or log scale
# ------------------------------------------------------------------------------
generate_heatmap <- function(data, title, scale_type = c("bidirectional", "log", "linear")) {
  scale_type <- match.arg(scale_type)
  
  # Base plot
  p <- ggplot(data, aes(x = proband1, y = proband2))
  
  if (scale_type == "bidirectional") {
  p <- p + 
    geom_tile(aes(fill = Relatedness)) +
    scale_y_discrete(limits = rev) +
    scale_fill_gradient2(
      low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
      trans = scales::pseudo_log_trans(sigma = 1e-4),  # Non-linear stretching
      na.value = "black",
      breaks = c(-1000, -100, -1, 0, 1, 100, 5000),
    )
}



 else if (scale_type == "log") {
    # Log scale for strictly positive data (with small offset)
    p <- p + 
      geom_tile(aes(fill = Relatedness + 1e-4)) +
      scale_y_discrete(limits = rev) +
      scale_fill_gradient(
        low = "white", high = "black", na.value = "black", trans = "log",
        breaks = 2^seq(-12, 0, by = 2),
        labels = format(2^seq(-12, 0, by = 2), digits = 3, scientific = FALSE)
      )
  } else {  # linear scale
    p <- p + 
      geom_tile(aes(fill = Relatedness)) +
      scale_y_discrete(limits = rev) +
      scale_fill_gradient(low = "white", high = "black", na.value = "black")
  }
  
  p + theme_minimal() +
    labs(title = title, x = "Individuals", y = "Individuals") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      legend.key.height = unit(0.4, "cm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}


# ------------------------------------------------------------------------------
# 5) Region color bars (manual color scheme)
# ------------------------------------------------------------------------------
generate_region_bar <- function(data) {
  # Define your custom palette
  region_colors <- c(
    "Chaleur Bay"   = "#E41A1C",
    "Batiscan"      = "#377EB8",
    "Chaudière"     = "#4DAF4A",
    "L'Assomption"  = "#984EA3",
    "Mistassini"    = "#FF7F00"
  )
  
  ggplot(
    distinct(data, proband1, proband_region1, depth1),
    aes(x = proband1, y = depth1, fill = proband_region1)
  ) +
    geom_col(color = NA) +
    scale_fill_manual(values = region_colors, name = "Region", na.value = "blue", drop = FALSE) +
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
# Process the GRM
res_grm <- process_relatedness(grm_file, meta_file)

# Process the PRM
res_prm <- process_relatedness(prm_file, meta_file)

# Use the PRM-based ranking to align both heatmaps
common_proband_ranking <- res_prm$ranking

# Extract the data frames
relatedness_grm <- res_grm$data
relatedness_prm <- res_prm$data

# Apply the common ranking
relatedness_grm <- apply_common_ranking(relatedness_grm, common_proband_ranking)
relatedness_prm <- apply_common_ranking(relatedness_prm, common_proband_ranking)

# Generate heatmaps:
# - GRM with linear scale
# - PRM with log scale
heatmap_grm <- generate_heatmap(
  relatedness_grm, "Estimated Genetic Relatedness Matrix (eGRM)", scale_type = "bidirectional"
)

heatmap_prm <- generate_heatmap(
  relatedness_prm, "Pedigree Relatedness Matrix (PRM)", scale_type = "log"
)

# Generate the region-color bars
region_colors_grm <- generate_region_bar(relatedness_grm)
region_colors_prm <- generate_region_bar(relatedness_prm)

# Combine heatmap and color bar for GRM
heatmap_grm_combined <- plot_grid(
  heatmap_grm + theme(plot.margin = margin(5, 0, 0, 5)),
  region_colors_grm + theme(plot.margin = margin(0, 0, 0, 5)),
  ncol = 1,
  rel_heights = c(1, 0.1),
  align = "v", 
  axis = "lr"
)

# Combine heatmap and color bar for PRM
heatmap_prm_combined <- plot_grid(
  heatmap_prm + theme(plot.margin = margin(5, 0, 0, 5)),
  region_colors_prm + theme(plot.margin = margin(0, 0, 0, 5)),
  ncol = 1,
  rel_heights = c(1, 0.1),
  align = "v", 
  axis = "lr"
)

# Combine both heatmaps, with PRM first and eGRM second
combined_plot <- plot_grid(
  heatmap_prm_combined, heatmap_grm_combined,  # Flipped order
  ncol = 2, align = "h", axis = "lr"
)

# Save the plot
ggsave(plot_file, plot = combined_plot, width = 15, height = 7, dpi = 300)
cat("Heatmaps with aligned color bars saved to:", plot_file, "\n")
