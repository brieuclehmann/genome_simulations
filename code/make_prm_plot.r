#!/usr/bin/env Rscript

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tools)
library(cowplot)
library(viridis)
library(ggdendro) # for converting hclust to ggplot-compatible data

set.seed(345)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the input file is provided
if (length(args) < 2) {
  stop("Usage: script.R input_file.csv metadata_file.csv")
}

# Path to the input files
data_file <- args[1]
meta_file <- args[2]

# Extract filename without path and extension
filename <- basename(data_file)
filename_no_ext <- file_path_sans_ext(filename)

# Determine output file base name and location
output_dir <- dirname(data_file)
plot1_file <- file.path(output_dir, paste0(filename_no_ext, "_boxplot.jpg"))
plot2_file <- file.path(output_dir, paste0(filename_no_ext, "_heatmap_with_bar_dendro.jpg"))

# ─────────────────────────────────────────────────────────────────────────────
# 1) Read the input files
# ─────────────────────────────────────────────────────────────────────────────
relatedness_matrix <- read_csv(data_file, col_names = TRUE)
metadata <- read_csv(meta_file, col_names = TRUE)

# Suppose 'metadata' now has columns:
#   proband, proband_region, average_depth
# e.g.
#   Batiscan,8.49557522123894,90
# Here 90 is the proband ID, "Batiscan" is the region, 8.49 is average depth

# Define an order for region factors if desired
region_order <- c("L'Assomption", "Batiscan", "Chaudière", "Mistassini", "Chaleur Bay")

# Convert proband_region to a factor with a specific order
metadata <- metadata %>%
  rename(depth = average_depth) %>%      # rename average_depth => depth
  mutate(
    proband_region = factor(proband_region, levels = region_order)
  ) %>%
  distinct(proband, proband_region, depth) 
# 'distinct()' ensures no duplicate rows if your CSV had multiples

# ─────────────────────────────────────────────────────────────────────────────
# 2) Pre-process the relatedness matrix
# ─────────────────────────────────────────────────────────────────────────────

# Add row names as a column using column names as IDs
relatedness_matrix <- relatedness_matrix %>%
  mutate(proband1 = colnames(relatedness_matrix))

# Convert your relatedness to a distance measure
# (the script does 1/(rel + 1e-8), or you might do (1 - rel) if that is your preferred metric)
mat <- as.matrix(relatedness_matrix[, -ncol(relatedness_matrix)]) # remove 'proband1' col
rownames(mat) <- relatedness_matrix$proband1
#dist_mat <- as.dist(1 / (mat + 1e-8)) # an inverted-similarity measure
dist_mat <- as.dist(1 - mat)

# Perform hierarchical clustering on rows (and columns)
hc <- hclust(dist_mat, method = "average")

# Save the hierarchical cluster order alone
clustered_proband <- rownames(mat)[hc$order]

# Build a region-based plus cluster-based ordering
region_map <- metadata %>%
  distinct(proband, proband_region) %>%
  filter(!is.na(proband_region))  # if any region is NA, decide how to handle

# Use match(...) to find each proband's position in the cluster order
region_priority <- region_map %>%
  mutate(cluster_order = match(proband, clustered_proband)) %>%
  # First sort by region, then within region by that cluster order
  arrange(proband_region, cluster_order)

# This combined ordering is region first, then cluster
proband_ranking <- region_priority$proband


# ─────────────────────────────────────────────────────────────────────────────
# 3) Reshape into a long format
# ─────────────────────────────────────────────────────────────────────────────
relatedness_long <- relatedness_matrix %>%
  pivot_longer(
    cols = -proband1,
    names_to = "proband2",
    values_to = "Relatedness"
  ) %>%
  mutate(
    # If you want to mask diagonal or something:
    Relatedness = ifelse(proband1 == proband2, NA, Relatedness),
    proband1 = as.numeric(proband1),
    proband2 = as.numeric(proband2)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 4) Join region/depth info for each proband
# ─────────────────────────────────────────────────────────────────────────────
# We'll do two separate joins: once for proband1, once for proband2
relatedness_long <- relatedness_long %>%
  left_join(
    metadata,
    by = c("proband1" = "proband")
  ) %>%
  rename(
    proband_region1 = proband_region,
    depth1 = depth
  ) %>%
  left_join(
    metadata,
    by = c("proband2" = "proband")
  ) %>%
  rename(
    proband_region2 = proband_region,
    depth2 = depth
  )

# If you also have a separate 'metadata' file with generation, relationship, etc.,
# you'd join that as well. For now let's assume generation/relationship are in this
# same or a separate structure. If you have them, do something like:
# relatedness_long <- relatedness_long %>%
#    left_join(other_table_of_relationships) 
# Then filter out or keep a "relatives" subset if generation is not NA.

# ─────────────────────────────────────────────────────────────────────────────
# 5) Apply cluster ordering
# ─────────────────────────────────────────────────────────────────────────────
relatedness_long$proband1 <- factor(relatedness_long$proband1, levels = proband_ranking)
relatedness_long$proband2 <- factor(relatedness_long$proband2, levels = proband_ranking)

# (OPTIONAL) If your metadata or prior step includes a 'generation' or 'relationship' column,
# then we can define:
relatives <- relatedness_long %>% filter(!is.na(Relatedness)) # or generation?
# If you have generation/relationship columns, do:
#   relatives <- relatedness_long %>% filter(!is.na(generation))

# This example does a simple boxplot if we had a 'relationship' column.
# If not, skip the boxplot by relationship.

## For demonstration, let's assume we do a simple distribution plot:
boxplot <- ggplot(relatives, aes(x = "", y = Relatedness)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Distribution of Relatedness",
    x = NULL,
    y = "Relatedness"
  )

ggsave(plot1_file, plot = boxplot, width = 8, height = 6, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Create a color bar for proband1 region using the actual depth
# ─────────────────────────────────────────────────────────────────────────────
# We'll just use (proband1, depth1) for the bar plot, coloring by proband_region1

region_colors <- ggplot(
  distinct(relatedness_long, proband1, proband_region1, depth1),
  aes(x = proband1, y = depth1, fill = proband_region1)
) +
  geom_col(color = NA) +
  scale_fill_viridis_d(option = "turbo", name = "Region", na.value = "blue") +
  scale_y_reverse() +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1))

# ─────────────────────────────────────────────────────────────────────────────
# 7) Generate the relatedness heatmap
# ─────────────────────────────────────────────────────────────────────────────
heatmap <- ggplot(relatedness_long, aes(x = proband1, y = proband2, fill = Relatedness + 1e-4)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white", high = "black",
    na.value = 'blue', trans = "log",
    breaks = c(0.001, 0.01, 0.1, 0.5)
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank()
  ) +
  labs(x = "Individuals", y = "Individuals")

# ─────────────────────────────────────────────────────────────────────────────
# 8) Build a title grob
# ─────────────────────────────────────────────────────────────────────────────
#title_grob <- ggdraw() + 
#  draw_label(
#    "Heatmap of Relatedness Matrix", 
#    fontface = 'bold', 
#    hjust = 0.5
#  )

# ─────────────────────────────────────────────────────────────────────────────
# 9) Create the dendrogram
# ─────────────────────────────────────────────────────────────────────────────
#dendro_data <- ggdendro::dendro_data(hc)

#dendrogram_plot <- ggplot(segment(dendro_data)) +
#  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#  scale_y_reverse() +  # you could try trans="log10" if you want
#  theme_classic() +
#  theme(
#    axis.title = element_blank(),
#    axis.text = element_blank(),
#    axis.ticks = element_blank(),
#    plot.margin = margin(0,0,0,0)
#  ) +
#  scale_x_continuous(
#    expand = c(0,0), 
#    limits = c(0.5, length(proband_ranking) + 0.5)
#  )

# ─────────────────────────────────────────────────────────────────────────────
# 10) Combine plots into one figure
# ─────────────────────────────────────────────────────────────────────────────
combined_plot <- plot_grid(
  #dendrogram_plot + theme(plot.margin = margin(0,0,0,5)),
  heatmap + theme(plot.margin = margin(5,0,0,5)),
  region_colors + theme(plot.margin = margin(0,0,0,5)),
  ncol = 1, 
  rel_heights = c(1, 0.1),
  align = "v", 
  axis = "lr"
)

#final_plot <- plot_grid(
#  title_grob,
#  combined_plot,
#  ncol = 1,
#  rel_heights = c(0.1, 1)
#)

# Save the combined plot
ggsave(plot2_file, plot = combined_plot, width = 8, height = 7, dpi = 300)

cat("Boxplot saved to:", plot1_file, "\n")
cat("Heatmap with region color bar and dendrogram saved to:", plot2_file, "\n")
