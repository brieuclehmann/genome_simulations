#!/usr/bin/env Rscript

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tools)
library(cowplot)
library(viridis)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the input file is provided
if (length(args) < 1) {
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
plot2_file <- file.path(output_dir, paste0(filename_no_ext, "_heatmap_with_bar.jpg"))

# Read the input files
relatedness_matrix <- read_csv(data_file, col_names = TRUE)
metadata <- read_csv(meta_file, col_names = TRUE)
region_order <-  c("L'Assomption", "Batiscan", "Chaudière", "Mistassini", "Chaleur Bay")

# Process metadata to extract proband regions
proband_region1 <- metadata %>%
  distinct(proband1, proband_region1) %>%
  rename(proband = proband1, proband_region = proband_region1)

proband_region2 <- metadata %>%
  distinct(proband2, proband_region2) %>%
  rename(proband = proband2, proband_region = proband_region2)

proband_regions <- bind_rows(proband_region1, proband_region2) %>%
  distinct(proband, proband_region)

# Add row names as a column using column names as IDs
relatedness_matrix <- relatedness_matrix %>%
  mutate(proband1 = colnames(relatedness_matrix))

proband_regions$proband_region <- factor(proband_regions$proband_region, levels = region_order)

# Transform data to long format
relatedness_long <- relatedness_matrix %>%
  pivot_longer(
    cols = -proband1,
    names_to = "proband2",
    values_to = "Relatedness"
  ) %>%
  mutate(Relatedness = ifelse(proband1 == proband2, NA, Relatedness),
         proband1 = as.numeric(proband1),
         proband2 = as.numeric(proband2)) %>%
  left_join(proband_regions, by = c("proband1" = "proband")) %>%
  left_join(proband_regions, by = c("proband2" = "proband"), suffix = c("1", "2")) %>%
  left_join(metadata) %>%
  mutate(
    proband_region1 = factor(proband_region1, levels =region_order),
    proband_region2 = factor(proband_region2, levels = region_order)
  )

# Filter relatives and set factor levels for relationship
relatives <- relatedness_long %>% filter(!is.na(generation))
rank <- relatives %>%
  distinct(generation, relationship) %>%
  arrange(generation) %>%
  pull(relationship)
relatives$relationship <- factor(relatives$relationship, levels = rank)

# Generate boxplot of relatedness by relationship
boxplot <- ggplot(relatives, aes(x = relationship, y = Relatedness)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme_bw() +
  labs(
    title = "Boxplot of Relatedness by Relationship",
    x = "Relationship",
    y = "Relatedness"
  )

# Save the boxplot
ggsave(plot1_file, plot = boxplot, width = 8, height = 6, dpi = 300)

# Create proband ranking for heatmap ordering
proband_ranking <- relatedness_long %>%
  arrange(proband_region1) %>%
  distinct(proband1) %>%
  pull(proband1)
relatedness_long$proband1 <- factor(relatedness_long$proband1, levels = proband_ranking)
relatedness_long$proband2 <- factor(relatedness_long$proband2, levels = proband_ranking)

# Create a dummy bar plot to indicate proband region
region_colors <- ggplot(relatedness_long,
                        aes(x = proband1, y = 1, fill = as.factor(proband_region1))) +
    geom_tile() +
    scale_fill_viridis_d(option = "turbo", name = "Region", na.value="blue") +
    theme_void() +
    theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
    guides(fill = guide_legend(nrow = 1))  # Force legend to one row


# Generate heatmap of relatedness
heatmap <- ggplot(relatedness_long, aes(x = proband1, y = proband2, fill = Relatedness + 1e-04)) +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Heatmap of Relatedness Matrix\nchr3 prm",
    x = "Individuals",
    y = "Individuals"
  )

# Combine the heatmap and color bar using cowplot
combined_plot <- plot_grid(
  heatmap + theme(plot.margin = margin(5, 0, 0, 0)),  
  region_colors + theme(plot.margin = margin(0, 0, 5, 0)),
  ncol = 1, rel_heights = c(1, 0.1), align = "v", axis = "lr"
)

# Save the combined plot
ggsave(plot2_file, plot = combined_plot, width = 8, height = 9, dpi = 300)

# Print the output file locations
cat("Boxplot saved to:", plot1_file, "\n")
cat("Heatmap with region color bar saved to:", plot2_file, "\n")
