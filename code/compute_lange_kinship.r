#!/usr/bin/env Rscript

# Computes genealogical (Lange) kinship and relationship labels up to k generations.
#
# Usage:
#   Rscript compute_lange_kinship.R pedigree.csv k [output.csv]
#
# Where:
#   pedigree.csv: A file with columns (ind, father, mother).
#   k: max number of generations to climb in BFS.
#   output.csv (optional): location to save the resulting table of:
#       proband1, proband2, kinship, relationship
#
# Example:
#   Rscript compute_lange_kinship.R my_pedigree.csv 4 results.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

###############################################################################
# 1) Enumerate all ancestor paths (BFS, up to k)
###############################################################################

enumerate_paths_up_to_k <- function(pedigree, focal_ind, k) {
  # BFS that records every path from focal_ind up to k edges (father/mother).
  #
  # Returns a data frame with columns:
  #   proband, ancestor, path_length, path_string
  #
  # path_string is mostly for debugging (e.g. "5->2->1").

  father_map <- setNames(pedigree$father, pedigree$ind)
  mother_map <- setNames(pedigree$mother, pedigree$ind)

  results <- data.frame(
    proband     = integer(),
    ancestor    = integer(),
    path_length = integer(),
    path_string = character(),
    stringsAsFactors = FALSE
  )

  # We'll keep a queue of lists: each has (current_id, length, path_str)
  queue <- list()
  queue[[1]] <- list(
    current_id = focal_ind,
    length = 0L,
    path_str = as.character(focal_ind)
  )
  q_idx <- 1

  while (q_idx <= length(queue)) {
    item <- queue[[q_idx]]
    q_idx <- q_idx + 1

    curr_id <- item$current_id
    curr_len <- item$length
    curr_path <- item$path_str

    # Record this path
    results <- rbind(
      results,
      data.frame(
        proband = focal_ind,
        ancestor = curr_id,
        path_length = curr_len,
        path_string = curr_path,
        stringsAsFactors = FALSE
      )
    )

    if (curr_len < k) {
      fa <- father_map[[as.character(curr_id)]]
      mo <- mother_map[[as.character(curr_id)]]
      for (parent_id in c(fa, mo)) {
        if (!is.na(parent_id) && parent_id != 0) {
          new_path_str <- paste0(curr_path, "->", parent_id)
          queue[[length(queue) + 1]] <- list(
            current_id = parent_id,
            length = curr_len + 1L,
            path_str = new_path_str
          )
        }
      }
    }
  }

  results
}

collect_paths_for_probands <- function(pedigree, probands, k) {
  # Enumerate BFS paths for each proband, combine into one data frame.
  do.call(rbind, lapply(probands, function(p) {
    enumerate_paths_up_to_k(pedigree, p, k)
  }))
}

###############################################################################
# 2) Compute Lange kinship by summing (1/2)^(path_len1 + path_len2)
###############################################################################

compute_lange_kinship_up_to_k <- function(pedigree, probands, k) {
  # 1) Gather all BFS paths
  all_paths <- collect_paths_for_probands(pedigree, probands, k)
  
  # We'll do a self-join on ancestor
  paths_l <- all_paths %>%
    rename(
      proband1 = proband,
      path_len1 = path_length,
      path_str1 = path_string
    )
  paths_r <- all_paths %>%
    rename(
      proband2 = proband,
      path_len2 = path_length,
      path_str2 = path_string
    )

  joined <- paths_l %>%
    inner_join(paths_r, by = "ancestor", relationship = "many-to-many") %>%
    filter(proband1 < proband2) %>%
    mutate(
      pair_id = paste0(proband1, "_", proband2),
      contrib = (1/2)^(path_len1 + path_len2)
    )

  # Sum up contributions
  kin_df <- joined %>%
    group_by(pair_id, proband1, proband2) %>%
    summarise(
      kinship = sum(contrib),
      .groups = "drop"
    ) %>%
    arrange(proband1, proband2)

  # Return both the summed kinship table and the joined path-level data
  list(kinship = kin_df, joined_paths = joined)
}

###############################################################################
# 3) Relationship labeling based on minimal path sum
###############################################################################

# We'll find the minimal path_length sum = (path_len1 + path_len2) for each pair.
# If sum=2 => siblings. Then we check father/mother in the pedigree for full vs. half.
# If sum=4 => "first_cousins"
# If sum=6 => "second_cousins"
# Otherwise => (sum/2 - 1)_degree_relatives (e.g. sum=8 => "3_degree_relatives")
# We skip half-cousins or more complicated cousin logic.

label_relationships <- function(joined_paths, pedigree) {
  joined_paths <- joined_paths %>%
    mutate(total_path_len = path_len1 + path_len2)

  # Minimal sum per pair
  min_path_summary <- joined_paths %>%
    group_by(proband1, proband2) %>%
    summarise(
      min_sum = min(total_path_len),
      .groups = "drop"
    )

  # Relationship name from sum
  classify_base_relationship <- function(s) {
    if (s == 2) {
      return("siblings")
    } else if (s == 4) {
      return("first_cousins")
    } else if (s == 6) {
      return("second_cousins")
    } else if (s > 2 && (s %% 2) == 0) {
      # For even sums >= 8 => (s/2 - 1)_degree_relatives
      g <- (s / 2) - 1
      return(paste0(g, "_degree_relatives"))
    } else {
      return("unclassified")
    }
  }

  # Sibling half/full check
  label_sibling_type <- function(p1, p2, ped) {
    # father, mother for each
    f1 <- ped$father[ped$ind == p1]
    m1 <- ped$mother[ped$ind == p1]
    f2 <- ped$father[ped$ind == p2]
    m2 <- ped$mother[ped$ind == p2]

    # if both father, mother match => full
    # if only one => half
    # else => "siblings"
    father_match <- (!is.na(f1) && !is.na(f2) && f1 == f2 && f1 != 0)
    mother_match <- (!is.na(m1) && !is.na(m2) && m1 == m2 && m1 != 0)

    if (father_match && mother_match) {
      return("full_siblings")
    } else if (father_match || mother_match) {
      return("half_siblings")
    } else {
      return("siblings")
    }
  }

  min_path_summary <- min_path_summary %>%
    mutate(
      base_relationship = sapply(min_sum, classify_base_relationship),
      final_relationship = base_relationship
    )

  # Upgrade siblings => full/half_siblings
  min_path_summary <- min_path_summary %>%
    rowwise() %>%
    mutate(
      final_relationship = if (base_relationship == "siblings") {
        label_sibling_type(proband1, proband2, pedigree)
      } else {
        base_relationship
      }
    ) %>%
    ungroup()

  min_path_summary
}

###############################################################################
# 4) Main script
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript compute_lange_kinship.R pedigree.csv k [output.csv]\n")
  quit(status = 1)
}

ped_file <- args[1]
k <- as.integer(args[2])
output_csv <- if (length(args) >= 3) args[3] else NULL

# Read pedigree
pedigree <- read.csv(ped_file, stringsAsFactors = FALSE)
pedigree <- pedigree %>%
  mutate(
    father = ifelse(father == 0, NA, father),
    mother = ifelse(mother == 0, NA, mother)
  )

# By default, define probands as those not listed as father/mother
all_parents <- c(pedigree$father, pedigree$mother)
probands <- setdiff(pedigree$ind, all_parents[!is.na(all_parents)])

# 1) Compute Lange kinship
res <- compute_lange_kinship_up_to_k(pedigree, probands, k)
kin_df <- res$kinship             # columns: pair_id, proband1, proband2, kinship
joined_paths <- res$joined_paths  # path-level detail

# 2) Label relationships from minimal path sums
rel_df <- label_relationships(joined_paths, pedigree)
# rel_df has columns: (proband1, proband2, min_sum, base_relationship, final_relationship)

# 3) Merge kinship + final relationship label
out_df <- kin_df %>%
  left_join(
    rel_df %>% select(proband1, proband2, final_relationship),
    by = c("proband1", "proband2")
  ) %>%
  select(proband1, proband2, kinship, relationship = final_relationship) %>%
  arrange(proband1, proband2)

cat(sprintf("Lange kinship among probands (up to %d generations) + relationship labels:\n", k))
print(out_df, n = 50)

# If an output CSV was provided, write it
if (!is.null(output_csv)) {
  write.csv(out_df, output_csv, row.names = FALSE)
  cat(sprintf("Wrote results to %s\n", output_csv))
}
