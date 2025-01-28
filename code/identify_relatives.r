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

# Recursive function to get ancestors up to k generations, now also tracking spouse information
get_ancestors <- function(pedigree, inds, k, current_generation = 0, lineage = "") {
  if (k < 0 || length(inds) == 0) {
    return(data.frame(
      ancestor = integer(0),
      spouse = integer(0),
      generation = integer(0),
      sex = integer(0),
      lineage = character(0),
      father_id = integer(0),
      mother_id = integer(0)
    ))
  }
  
  # Parse input inds (may contain ancestorID_spouseID pairs)
  parsed <- lapply(inds, function(x) {
    if (grepl("_", x, fixed = TRUE)) {
      parts <- strsplit(x, "_", fixed = TRUE)[[1]]
      ancestor_id <- as.integer(parts[1])
      spouse_id <- ifelse(parts[2] == "" || parts[2] == "NA", NA, as.integer(parts[2]))
      list(ancestor = ancestor_id, spouse = spouse_id)
    } else {
      list(ancestor = as.integer(x), spouse = NA)
    }
  })
  
  parsed_df <- do.call(rbind, lapply(parsed, as.data.frame, stringsAsFactors = FALSE))
  
  # Retrieve information for current ancestors from pedigree
  current_info <- pedigree %>% filter(ind %in% parsed_df$ancestor) %>%
    select(ind, sex, father, mother)
  
  # Add spouse information from parsed_df
  current_ancestors <- current_info %>%
    left_join(parsed_df, by = c("ind" = "ancestor")) %>%
    mutate(
      generation = current_generation,
      lineage = lineage,
      father_id = ifelse(father == 0, NA, father),
      mother_id = ifelse(mother == 0, NA, mother)
    ) %>%
    select(ancestor = ind, spouse, generation, sex, lineage, father_id, mother_id)
  
  # Generate next generation pairs with ancestor and spouse information
  next_pairs <- list()
  for (i in seq_len(nrow(current_ancestors))) {
    row <- current_ancestors[i, ]
    fa <- row$father_id
    mo <- row$mother_id
    
    # Handle NA parents by substituting "NA" in the pair string
    fa_str <- ifelse(is.na(fa), "NA", as.character(fa))
    mo_str <- ifelse(is.na(mo), "NA", as.character(mo))
    
    if (!(is.na(fa) && is.na(mo))) {
      next_pairs[[length(next_pairs) + 1]] <- paste0(mo_str, "_", fa_str)
      next_pairs[[length(next_pairs) + 1]] <- paste0(fa_str, "_", mo_str)
    }
  }
  
  next_inds <- unique(unlist(next_pairs))
  
  # Recurse upwards
  upper_founders <- get_ancestors(
    pedigree,
    next_inds,
    k - 1,
    current_generation + 1,
    lineage = lineage
  )
  
  all_ancestors <- bind_rows(current_ancestors, upper_founders) %>%
    distinct()
  
  return(all_ancestors)
}

# Function to get all ancestors for each proband up to k generations, parallelized
get_all_ancestors_furrr <- function(pedigree, list_of_probands, k) {
  plan(multisession, workers = max(1, future::availableCores() - 1))

  all_ancestors <- future_map_dfr(list_of_probands, function(proband) {
    ancestors <- get_ancestors(pedigree, proband, k)
    ancestors$proband <- proband
    return(ancestors)
  })

  all_ancestors <- all_ancestors %>% filter(!is.na(ancestor))

  return(all_ancestors)
}

# Function to identify common ancestors among probands
identify_common_ancestors <- function(all_ancestors) {
  common_ancestors <- all_ancestors %>%
    group_by(ancestor, generation) %>%
    filter(n_distinct(proband) > 1) %>%
    ungroup()

  return(common_ancestors)
}

# Function to group relatives based on shared ancestors and their spouses
group_relatives <- function(common_ancestors) {
  relatives <- common_ancestors %>%
    mutate(relationship = case_when(
      generation == 1 ~ "siblings",
      generation == 2 ~ "first_cousins",
      generation == 3 ~ "second_cousins",
      TRUE ~ paste0(generation - 1, "_degree_relatives")
    ))

  # Create all pair combinations of probands for each common ancestor
  relative_pairs <- relatives %>%
    select(proband, ancestor, generation, sex, lineage, father_id, mother_id, spouse, relationship) %>%
    group_by(ancestor, generation, relationship) %>%
    summarise(
      pairs = list(as.data.frame(t(combn(proband, 2)))),
      .groups = "drop"
    ) %>%
    unnest_longer(pairs) %>%
    unnest_wider(pairs) %>%
    rename(proband1 = V1, proband2 = V2)

  symmetric_pairs <- relative_pairs %>%
    mutate(temp = proband1, proband1 = proband2, proband2 = temp) %>%
    select(-temp)

  all_pairs <- bind_rows(relative_pairs, symmetric_pairs) %>% distinct()

  # Join back to get lineage, spouse info for both probands
  all_pairs <- all_pairs %>%
    left_join(
      relatives %>% select(proband, ancestor, generation, relationship, lineage, spouse),
      by = c("proband1" = "proband", "ancestor", "generation", "relationship")
    ) %>%
    rename(lineage1 = lineage, spouse1 = spouse) %>%
    left_join(
      relatives %>% select(proband, ancestor, generation, relationship, lineage, spouse),
      by = c("proband2" = "proband", "ancestor", "generation", "relationship")
    ) %>%
    rename(lineage2 = lineage, spouse2 = spouse)

  return(all_pairs)
}

# Main function
identify_related_probands <- function(pedigree, list_of_probands, k) {
  all_ancestors <- get_all_ancestors_furrr(pedigree, list_of_probands, k)
  common_ancestors <- identify_common_ancestors(all_ancestors)
  relatives <- group_relatives(common_ancestors)
  return(relatives)
}

###############################################################
# MAIN SCRIPT
###############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript identify_relatives.R input_pedigree.csv k [output_file.csv]\n")
  quit(status = 1)
}

input_file <- args[1]
k_value <- as.numeric(args[2])
output_file <- ifelse(length(args) == 3, args[3], "")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

if (is.na(k_value) || k_value < 1) {
  stop("k must be a positive integer")
}

pedigree <- data.table::fread(input_file) %>%
  mutate(father = ifelse(father == 0, NA, father),
         mother = ifelse(mother == 0, NA, mother)) %>%
  filter(mother != 0 | is.na(mother), father != 0 | is.na(father), ind != 0)

list_of_probands <- pedigree %>% 
  filter(!ind %in% c(pedigree$mother, pedigree$father)) %>% 
  pull(ind) 

relatives <- identify_related_probands(pedigree, list_of_probands, k_value)

# Keep only the closest relationship for each pair
relatives <- relatives %>%
  group_by(proband1, proband2) %>%
  filter(generation == min(generation)) %>%
  ungroup() %>%
  distinct()

# Determine relationship type from (ancestor, spouse)
# If both probands share exactly the same (ancestor, spouse), we have full.
# If they share the ancestor but differ in spouse, we have half.
# If multiple distinct pairs appear, other.
# Ensure that spouse1 and spouse2 are joined before the mutate step
relatives_with_type <- relatives %>%
  left_join(
    relatives %>% select(proband, ancestor, generation, relationship, spouse) %>%
      rename(spouse1 = spouse),
    by = c("proband1" = "proband", "ancestor", "generation", "relationship")
  ) %>%
  left_join(
    relatives %>% select(proband, ancestor, generation, relationship, spouse) %>%
      rename(spouse2 = spouse),
    by = c("proband2" = "proband", "ancestor", "generation", "relationship")
  ) %>%
  group_by(ancestor, generation, relationship, proband1, proband2) %>%
  summarise(
    distinct_pairs = n_distinct(paste0(
      ifelse(is.na(spouse1), "NA", spouse1), "_",
      ifelse(is.na(spouse2), "NA", spouse2)
    )),
    both_parents_present = all(!is.na(father_id1), !is.na(mother_id1), !is.na(father_id2), !is.na(mother_id2)),
    .groups = "drop"
  ) %>%
  mutate(relationship_type = case_when(
    distinct_pairs == 1 & both_parents_present ~ "full",
    distinct_pairs == 1 ~ "half",
    distinct_pairs > 1  ~ "other"
  )) %>%
  left_join(relatives, by = c("ancestor", "generation", "relationship", "proband1", "proband2"))

relatives_with_type <- relatives_with_type %>% arrange(generation, ancestor)


cat("Counts of relationship types:\n")
print(table(relatives_with_type$relationship_type))

if (output_file != "") {
  write.csv(relatives_with_type, output_file, row.names = FALSE)
  cat("Results saved to:", output_file, "\n")
} else {
  print(relatives_with_type)
}
