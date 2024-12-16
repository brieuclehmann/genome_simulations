#!/usr/bin/env Rscript

# Load required packages
library(dplyr)
library(tidyr)
library(data.table)

# Usage:
# Rscript find_founder_depths.R input_pedigree.csv [output_file.csv]
#
# If output_file.csv is provided, the results are written to that file and only the header is printed to stdout.
# If no output_file.csv is provided, results are printed to stdout.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript find_founder_depths.R input_pedigree.csv [output_file.csv]\n")
  quit(status = 1)
}

input_file <- args[1]
output_file <- if (length(args) == 2) args[2] else NULL

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

# Read and preprocess the pedigree data
# We assume a standard structure with columns: ind, father, mother
pedigree <- fread(input_file)

# Ensure that missing parents are represented as NA
pedigree <- pedigree %>%
  mutate(
    father = ifelse(father == 0, NA, father),
    mother = ifelse(mother == 0, NA, mother)
  )

# Identify probands (individuals without recorded offspring)
# These are individuals that do not appear as father/mother of anyone.
list_of_probands <- pedigree %>%
  filter(!ind %in% c(pedigree$father, pedigree$mother)) %>%
  pull(ind)

# Recursive function to get founder depths
get_founder_depths <- function(pedigree, inds, current_generation = 0) {
  if (length(inds) == 0) {
    # No individuals provided, return empty
    return(data.frame(founder = integer(0), depth = integer(0)))
  }
  
  # For each individual, retrieve their parents
  parents <- pedigree %>% filter(ind %in% inds) %>%
    select(ind, father, mother)
  
  # Check which are founders: those with no known parents
  founders <- parents %>%
    filter(is.na(father) & is.na(mother)) %>%
    select(founder = ind)
  
  # If we found some founders, record their depths
  founder_depths <- data.frame()
  if (nrow(founders) > 0) {
    founder_depths <- founders %>%
      mutate(depth = current_generation)
  }
  
  # For non-founders, recursively climb up
  non_founders <- parents %>%
    filter(!(is.na(father) & is.na(mother)))
  
  if (nrow(non_founders) > 0) {
    # Collect all parent IDs that are not missing
    next_inds <- c(na.omit(non_founders$father), na.omit(non_founders$mother))
    # Recurse one generation up
    upper_founders <- get_founder_depths(pedigree, next_inds, current_generation = current_generation + 1)
    # Combine with current founders
    founder_depths <- bind_rows(founder_depths, upper_founders) %>% distinct()
  }

  return(founder_depths)
}

# For each proband, find the founder depths and compute the average
results <- lapply(list_of_probands, function(p) {
  df <- get_founder_depths(pedigree, p, current_generation = 0)
  if (nrow(df) == 0) {
    # No founders found, unusual, but handle gracefully
    return(data.frame(
      proband = p,
      average_depth = NA,
      founder_ids = I(list(integer(0))),
      founder_depths = I(list(integer(0)))
    ))
  } else {
    avg_depth <- mean(df$depth)
    return(data.frame(
      proband = p,
      average_depth = avg_depth,
      founder_ids = I(list(df$founder)),
      founder_depths = I(list(df$depth))
    ))
  }
}) %>% bind_rows()

# Print only the header if output file is specified, otherwise print the full results
if (!is.null(output_file)) {
  # Write results to file
  write.csv(results, output_file, row.names = FALSE)
  
  # Print only the header (column names) to stdout
  cat(paste(names(results), collapse = ","), "\n")
} else {
  # Print full results to stdout
  print(results)
}
