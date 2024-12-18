#!/usr/bin/env Rscript

# Load required packages
library(dplyr)
library(tidyr)
library(furrr)
library(purrr)
library(data.table)
library(future)

###############################################################
# FUNCTIONS
###############################################################

# Recursive function to get ancestors up to k generations, including sex and lineage
get_ancestors <- function(pedigree, inds, k, current_generation = 0, lineage = "") {
  if (k < 0 || length(inds) == 0) {
    # Base case: Return an empty data frame if no individuals or generations left
    return(data.frame(ancestor = integer(0), generation = integer(0), sex = integer(0), lineage = character(0)))
  } else {
    # Current ancestors at this generation with lineage information
    current_ancestors <- data.frame(ancestor = inds, generation = current_generation, lineage = lineage) %>%
      left_join(pedigree %>% select(ind, sex), by = c("ancestor" = "ind"))
    
    # Retrieve parents of the current individuals
    parents <- pedigree %>%
      filter(ind %in% inds) %>%
      select(father, mother)
    
    # Recursively get ancestors from father's side, appending 'p' to lineage
    ancestors_father <- get_ancestors(
      pedigree,
      na.omit(parents$father),
      k - 1,
      current_generation + 1,
      lineage = paste(lineage, "p", sep = ">")
    )
    
    # Recursively get ancestors from mother's side, appending 'm' to lineage
    ancestors_mother <- get_ancestors(
      pedigree,
      na.omit(parents$mother),
      k - 1,
      current_generation + 1,
      lineage = paste(lineage, "m", sep = ">")
    )
    
    # Combine all ancestors and remove duplicates
    all_ancestors <- bind_rows(current_ancestors, ancestors_father, ancestors_mother) %>%
      distinct()
    
    return(all_ancestors)
  }
}

# Function to get all ancestors for each proband up to k generations, parallelized with furrr
get_all_ancestors_furrr <- function(pedigree, list_of_probands, k) {
  # Plan for parallel processing using available cores
  plan(multisession, workers = max(1, future::availableCores() - 1))
  
  # Use future_map_dfr for parallel execution over probands
  all_ancestors <- future_map_dfr(list_of_probands, function(proband) {
    ancestors <- get_ancestors(pedigree, proband, k)
    ancestors$proband <- proband
    return(ancestors)
  })
  
  # Remove missing ancestors
  all_ancestors <- all_ancestors %>% filter(!is.na(ancestor))
  
  return(all_ancestors)
}

# Function to identify common ancestors among probands
identify_common_ancestors <- function(all_ancestors) {
  # Identify ancestors shared by more than one proband
  common_ancestors <- all_ancestors %>%
    group_by(ancestor, generation) %>%
    filter(n_distinct(proband) > 1) %>%
    ungroup()
  
  return(common_ancestors)
}

# Function to group relatives based on shared ancestors
group_relatives <- function(common_ancestors) {
  # Define relationships based on the generation of the shared ancestor
  relatives <- common_ancestors %>%
    mutate(relationship = case_when(
      generation == 1 ~ "siblings",
      generation == 2 ~ "first_cousins",
      generation == 3 ~ "second_cousins",
      TRUE ~ paste0(generation - 1, "_degree_relatives")
    ))
  
  # Create all pair combinations of probands for each common ancestor
  relative_pairs <- relatives %>%
    select(proband, ancestor, generation, sex, lineage, relationship) %>%
    group_by(ancestor, generation, relationship) %>%
    summarise(
      pairs = list(as.data.frame(t(combn(proband, 2)))),
      .groups = "drop"
    ) %>%
    unnest_longer(pairs) %>%
    unnest_wider(pairs) %>%
    rename(proband1 = V1, proband2 = V2)
  
  # Add symmetry by swapping probands and duplicating rows
  symmetric_pairs <- relative_pairs %>%
    mutate(temp = proband1, proband1 = proband2, proband2 = temp) %>%
    select(-temp)
  
  # Combine original and symmetric pairs and remove duplicates
  all_pairs <- bind_rows(relative_pairs, symmetric_pairs) %>%
    distinct()
  
  # Add lineage for proband1 only (optional, may not be critical)
  all_pairs <- all_pairs %>%
    left_join(
      relatives %>% select(proband, ancestor, lineage) %>% rename(lineage1 = lineage),
      by = c("proband1" = "proband", "ancestor" = "ancestor")
    )
  
  return(all_pairs)
}

# Main function to identify and classify relative types among probands
identify_related_probands <- function(pedigree, list_of_probands, k) {
  # Step 1: Get all ancestors up to generation k
  all_ancestors <- get_all_ancestors_furrr(pedigree, list_of_probands, k)
  
  # Step 2: Identify common ancestors among probands
  common_ancestors <- identify_common_ancestors(all_ancestors)
  
  # Step 3: Group relatives based on shared ancestors and their generation
  relatives <- group_relatives(common_ancestors)
  
  # Return the result
  return(relatives)
}

###############################################################
# MAIN SCRIPT
###############################################################

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript identify_relatives.R input_pedigree.csv k [output_file.csv]\n")
  quit(status = 1)
}

input_file <- args[1]
k_value <- as.numeric(args[2])
output_file <- ifelse(length(args) == 3, args[3], "")

# Sanity checks
if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

if (is.na(k_value) || k_value < 1) {
  stop("k must be a positive integer")
}

# Read and preprocess the pedigree data
pedigree <- data.table::fread(input_file) %>% 
  filter(mother != 0, father != 0, ind != 0)

# Identify probands (individuals without recorded offspring in the pedigree)
list_of_probands <- pedigree %>% 
  filter(!ind %in% c(pedigree$mother, pedigree$father)) %>% 
  pull(ind) 

# Identify related probands up to k generations
relatives <- identify_related_probands(pedigree, list_of_probands, k_value) 

# Process relatives to keep only the closest relationship for each pair
relatives <- relatives %>%
  group_by(proband1, proband2) %>%
  filter(generation == min(generation)) %>%
  ungroup() %>%
  distinct()

# Arrange relatives by generation and ancestor
relatives <- relatives %>% arrange(generation, ancestor)

# Output results
if (output_file != "") {
  write.csv(relatives, output_file, row.names = FALSE)
  cat("Results saved to:", output_file, "\n")
} else {
  print(relatives)
}
