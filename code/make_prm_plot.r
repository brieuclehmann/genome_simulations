#!/usr/bin/env Rscript

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tools)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the input file is provided
if(length(args) < 1) {
  stop("Usage: script.R input_file.csv")
}

# Path to the CSV file
data_file <- args[1]
# Extract filename without path and extension
filename <- basename(data_file)
filename_no_ext <- file_path_sans_ext(filename)

# Determine output file name and location
if (length(args) == 2) {
  output_file <- args[2]  # Use specified output file
} else {
  output_file <- file.path(dirname(data_file), paste0(filename_no_ext, ".jpg"))
}

# Read the data (assuming headers)
relatedness_matrix <- read_csv(data_file, col_names = TRUE)

# Get the number of rows and columns
n_rows <- nrow(relatedness_matrix)
n_cols <- ncol(relatedness_matrix)

# Generate row and column names
row_names <- paste0("Ind", 1:n_rows)
col_names <- paste0("Ind", 1:n_cols)

# Assign column names
colnames(relatedness_matrix) <- col_names

# Transform data to long format and remove diagonal values
relatedness_long <- relatedness_matrix %>%
  mutate(Row = row_names) %>%
  pivot_longer(
    cols = starts_with("Ind"),
    names_to = "Column",
    values_to = "Value"
  ) %>% 
  mutate(Value = ifelse(Row == Column, NA, Value))

# Create plot title by replacing underscores with spaces
plot_title <- gsub("_", " ", filename_no_ext)
plot_title <- paste("Heatmap of Relatedness Matrix\n", plot_title)

# Create the heatmap
p <- ggplot(relatedness_long, aes(x = Column, y = Row, fill = -log(Value + 1e-04))) +
  geom_tile() +
  scale_fill_gradient(high = "white", low = "black", na.value = 'blue') +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    # Hide x-axis labels
    axis.text.y = element_blank(),    # Hide y-axis labels
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = plot_title,
    x = "Individuals",
    y = "Individuals"
  )

# Save the plot to a JPEG file
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)

# Print the output file location
cat("Heatmap saved to:", output_file, "\n")
