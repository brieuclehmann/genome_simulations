#!/usr/bin/env Rscript

# pca_plot_depth_scaled.R
# - Reads PCA data from branch and pedigree files
# - Reads metadata file for region and depth
# - Joins metadata with PCA data using `noid` and verifies region consistency
# - Generates PCA plots (PC1-2, PC3-4, PC5-6) with point size scaled inversely to depth
# - Saves the final combined figure
# - BELOW THIS, we also generate *heatmaps* with the requested modifications

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  # If you need the below packages for the heatmap portion, ensure they're installed:
  # library(cowplot)
  # library(viridis)
  # library(tools)
  # library(tibble)
  # library(scales)
})

# ─────────────────────────────────────────────────────────────────────────────
# 1) File Paths
# ─────────────────────────────────────────────────────────────────────────────

BasePath         <- "~/Documents/genome_simulations_tsrelatedness/misc/geo_map/"
pca_branch_csv   <- paste0(BasePath, "balsac_branch_pca_noshallow.csv")
pca_pedigree_csv <- paste0(BasePath, "balsac_pedigree_pca_noshallow.csv")
meta_csv         <- paste0(BasePath, "balsac_proband_meta_noid.csv")
output_file      <- paste0(BasePath, "pca_depth_scaled_noshallow.jpg")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Define region color palette
# ─────────────────────────────────────────────────────────────────────────────

region_colors <- c(
  "Chaleur Bay"   = "#E41A1C",
  "Batiscan"      = "#377EB8",
  "Chaudière"     = "#4DAF4A",
  "L’Assomption"  = "#984EA3",
  "Mistassini"    = "#FF7F00"
)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Load and Process Data
# ─────────────────────────────────────────────────────────────────────────────

# Function to load and process PCA data
load_pca_data <- function(file_path) {
  fread(file_path, header = TRUE) %>%
    rename(
      PC1 = `0`, PC2 = `1`, PC3 = `2`, PC4 = `3`,
      PC5 = `4`, PC6 = `5`, region_pca = proband_region, ID = noid
    ) %>%
    select(ID, PC1, PC2, PC3, PC4, PC5, PC6, region_pca) %>%
    mutate(region_pca = str_replace_all(region_pca, "'", "’"))
}

# Load PCA files
pca_branch_df   <- load_pca_data(pca_branch_csv)
pca_pedigree_df <- load_pca_data(pca_pedigree_csv)

# Load metadata file
meta_df <- fread(meta_csv, header = TRUE) %>%
  rename(region_meta = proband_region, depth = average_depth, ID = noid) %>%
  select(ID, region_meta, depth) %>%
  mutate(region_meta = str_replace_all(region_meta, "'", "’"))

# Print count of individuals with depth = 0, 1, or 2
depth_counts <- meta_df %>%
  filter(depth %in% c(0, 1, 2)) %>%
  group_by(depth) %>%
  summarise(count = n(), .groups = "drop")

cat("\n🔍 Depth Value Counts:\n")
print(depth_counts)

print(head(meta_df %>% arrange(depth) %>% as.data.frame(), 30))

# Merge PCA data with metadata
pca_branch_df   <- left_join(pca_branch_df, meta_df, by = "ID")
pca_pedigree_df <- left_join(pca_pedigree_df, meta_df, by = "ID")

# Check for mismatches
check_mismatches <- function(df, label) {
  mismatches <- sum(df$region_pca != df$region_meta, na.rm = TRUE)
  if (mismatches > 0) {
    cat("\n⚠️ Mismatches detected in", label, ":", mismatches, "\n")
    print(head(df[df$region_pca != df$region_meta, ], 5))
    stop("Region mismatch detected.")
  } else {
    cat("\n✅ No mismatches detected in", label, "\n")
  }
}

check_mismatches(pca_branch_df, "Branch PCA")
check_mismatches(pca_pedigree_df, "Pedigree PCA")

# Inverse scaling of depth for point size
epsilon <- 0.0001
pca_branch_df   <- pca_branch_df   %>% mutate(size = 1 / (depth + epsilon))
pca_pedigree_df <- pca_pedigree_df %>% mutate(size = 1 / (depth + epsilon))

print(head(pca_branch_df %>% arrange(-size) %>% 
            select(ID, PC1, region_pca, depth, size)%>% as.data.frame(), 30))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Generate PCA Plots
# ─────────────────────────────────────────────────────────────────────────────

create_pca_plot <- function(df, x_var, y_var, x_label, y_label) {
  ggplot(df %>% arrange(size), aes_string(x = x_var, y = y_var, 
                                          color = "region_pca", size = "size")) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = region_colors) +
    scale_size_area(max_size = 5) +
    labs(x = x_label, y = y_label) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14)
    )
}

# Create PCA plots
pca_plot_12 <- create_pca_plot(pca_branch_df, "PC1", "PC2", "PC1", "PC2")
pca_plot_34 <- create_pca_plot(pca_branch_df, "PC3", "PC4", "PC3", "PC4")
pca_plot_56 <- create_pca_plot(pca_branch_df, "PC5", "PC6", "PC5", "PC6")

pca_pedigree_plot_12 <- create_pca_plot(pca_pedigree_df, "PC1", "PC2", "PC1", "PC2")
pca_pedigree_plot_34 <- create_pca_plot(pca_pedigree_df, "PC3", "PC4", "PC3", "PC4")
pca_pedigree_plot_56 <- create_pca_plot(pca_pedigree_df, "PC5", "PC6", "PC5", "PC6")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Combine Plots (PCA) and Save
# ─────────────────────────────────────────────────────────────────────────────

combined_plot <- (pca_plot_12 | pca_plot_34 | pca_plot_56) /
                 (pca_pedigree_plot_12 | pca_pedigree_plot_34 | pca_pedigree_plot_56) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

# Save final combined figure with adjusted dimensions
ggsave(filename = output_file, plot = combined_plot, width = 15, height = 10, dpi = 300)
cat("\n✅ PCA depth-scaled plot saved to:", output_file, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 6) Heatmap Changes (Minimal Edits)
# ─────────────────────────────────────────────────────────────────────────────

# We assume you have the following or similar heatmap code somewhere else.
# This snippet shows how to:
#  (1) flip the y-axis
#  (2) remove axis text and reduce the margin
#  (3) use integer-like breaks for the PRM color scale on log scale

library(scales)  # for breaks/labels if needed

generate_heatmap <- function(data, title, scale_type = c("linear","log","asinh")) {
  scale_type <- match.arg(scale_type)
  
  # Start with flipping y-axis via scale_y_reverse()
  p <- ggplot(data, aes(x = proband1, y = proband2)) + 
    scale_y_reverse()
  
  if (scale_type == "linear") {
    p <- p + 
      geom_tile(aes(fill = Relatedness)) +
      scale_fill_gradient(low = "white", high = "black", na.value = "blue")
    
  } else if (scale_type == "log") {
    # Use nicer breakpoints for a log scale:
    p <- p + 
      geom_tile(aes(fill = Relatedness + 1e-4)) +
      scale_fill_gradient(
        low = "white", high = "black", na.value = "blue",
        trans = "log",
        # Example of integer or "nice" breaks for your range:
        breaks = c(1e-3, 1e-2, 1e-1, 1),
        labels = c("0.001", "0.01", "0.1", "1")
      )
    
  } else if (scale_type == "asinh") {
    asinh_trans <- trans_new("asinh", transform = asinh, inverse = sinh)
    p <- p + 
      geom_tile(aes(fill = Relatedness)) +
      scale_fill_gradient(
        low = "white", high = "black", na.value = "blue",
        trans = asinh_trans
      )
  }
  
  # Remove axis text, axis titles, and reduce margins
  p + 
    theme_minimal() +
    labs(title = title) +
    theme(
      # Remove axis texts and titles:
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks  = element_blank(),
      # Slightly reduce margins:
      plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(hjust = 0.5)
    )
}

