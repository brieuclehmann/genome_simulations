#!/usr/bin/env Rscript

# This script reads a pedigree CSV with columns:
#   ind, father, mother
# (with 0 or NA for unknown). It computes genealogical kinship
# (Lange definition) for all proband pairs, restricted to k generations.
#
# Usage example:
#   Rscript compute_lange_kinship.R pedigree.csv 3
#   # => enumerates all paths up to 3 generations, sums path-based kinship.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

###############################################################################
# 1) Enumerate all ancestor paths (up to k generations) for each proband
###############################################################################

enumerate_paths_up_to_k <- function(pedigree, focal_ind, k) {
  # We do a BFS/DFS approach that keeps track of every distinct path from
  # focal_ind up to an ancestor within k generations.
  #
  # INPUT:
  #   - pedigree: data.frame with (ind, father, mother), NA for unknown.
  #   - focal_ind: single individual ID (integer)
  #   - k: max number of parent-child steps to climb
  #
  # OUTPUT: data.frame with columns:
  #   proband      => focal_ind
  #   ancestor     => ID of the ancestor
  #   path_length  => how many "parent-child" links in the path
  #   path_string  => string describing the path (optional), e.g. "child->father->..."
  #
  #   Each row is a unique path from focal_ind to 'ancestor'.

  # Create a quick lookup for father,mother by ID
  # This speeds up parent lookups
  father_map <- setNames(pedigree$father, pedigree$ind)
  mother_map <- setNames(pedigree$mother, pedigree$ind)

  # We'll store results in a data frame
  results <- data.frame(
    proband = integer(),
    ancestor = integer(),
    path_length = integer(),
    path_string = character(),
    stringsAsFactors = FALSE
  )

  # We'll keep a queue of paths. Each path is a list with:
  #   current_id: current node
  #   length: number of edges from focal_ind
  #   path_str: a textual representation (optional)
  queue <- list()
  queue[[1]] <- list(
    current_id = focal_ind,
    length = 0L,
    path_str = as.character(focal_ind)
  )
  q_idx <- 1

  visited_paths <- list()  # optional if you want to avoid exact path duplicates

  # BFS/DFS loop
  while (q_idx <= length(queue)) {
    item <- queue[[q_idx]]
    q_idx <- q_idx + 1

    curr_id <- item$current_id
    curr_len <- item$length
    curr_path <- item$path_str

    # Record the path in results
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

    # If we haven't reached k steps, expand father and mother
    if (curr_len < k) {
      fa <- father_map[[as.character(curr_id)]]
      mo <- mother_map[[as.character(curr_id)]]

      for (parent_id in c(fa, mo)) {
        if (!is.na(parent_id) && parent_id != 0) {
          new_path_str <- paste0(curr_path, "->", parent_id)
          # We could check if we've already enqueued this exact path,
          # but to allow distinct routes, we won't block duplicates here.
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
  # For each proband, enumerate all paths up to k generations.
  # Return a single data.frame with columns: proband, ancestor, path_length, path_string
  all_paths <- lapply(probands, function(p) {
    df_p <- enumerate_paths_up_to_k(pedigree, p, k)
    df_p
  })
  do.call(rbind, all_paths)
}


###############################################################################
# 2) Compute Lange kinship
###############################################################################
#
# We define the genealogical (Lange) kinship coefficient as the sum over
# all distinct paths from i -> A and j -> A of (1/2)^(path_length_i + path_length_j),
# for any shared ancestor A. Here "distinct" means distinct BFS paths in each proband's data.
#
# Note that if a proband can reach the same ancestor multiple ways, each path is counted.

compute_lange_kinship_up_to_k <- function(pedigree, probands, k) {
  # 1) Enumerate all paths for each proband
  all_paths <- collect_paths_for_probands(pedigree, probands, k)
  # The result has:
  #    proband, ancestor, path_length, path_string

  # 2) We'll do a "self join" on (ancestor) to find pairs of proband paths that share that ancestor
  #    Then sum over (1/2)^(path_length_i + path_length_j).
  #    We want all pairs (p1, p2) with p1 < p2 to avoid double counting.

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
    inner_join(paths_r, by = "ancestor") %>%
    filter(proband1 < proband2)  # ensure proband1 < proband2

  # Now each row in 'joined' corresponds to a distinct path from proband1 -> ancestor
  # and a distinct path from proband2 -> ancestor. The combined path length is sum(path_len1, path_len2).
  # The contribution to kinship = (1/2)^(path_len1 + path_len2).

  joined <- joined %>%
    mutate(
      pair_id = paste0(proband1, "_", proband2),
      contrib = (1/2)^(path_len1 + path_len2)
    )

  # 3) Summarize by pair_id to get total Lange kinship
  kin_df <- joined %>%
    group_by(pair_id, proband1, proband2) %>%
    summarise(
      kinship = sum(contrib),
      .groups = "drop"
    ) %>%
    arrange(proband1, proband2)

  # If you want a complete matrix including p1 == p2 or pairs that do not share ancestors,
  # you can expand the grid. But typically we just output the pairs that have a > 0 kinship.
  kin_df
}


###############################################################################
# 3) Main driver
###############################################################################

# Minimal CLI usage:
#   Rscript compute_lange_kinship.R pedigree.csv 3
#   => prints kinship among probands up to 3 generations

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript compute_lange_kinship.R pedigree.csv k\n")
  quit(status = 1)
}

ped_file <- args[1]
k <- as.integer(args[2])

# Read the pedigree
pedigree <- read.csv(ped_file, stringsAsFactors = FALSE)

# If father or mother are 0 => treat as NA
pedigree <- pedigree %>%
  mutate(
    father = ifelse(father == 0, NA, father),
    mother = ifelse(mother == 0, NA, mother)
  )

# By default, define probands as those who do not appear as father or mother
all_parents <- c(pedigree$father, pedigree$mother)
probands <- setdiff(pedigree$ind, all_parents[!is.na(all_parents)])

# Compute genealogical kinship
kin_df <- compute_lange_kinship_up_to_k(pedigree, probands, k)

cat(sprintf("Genealogical (Lange) kinship among probands, up to %d generations:\n", k))
print(kin_df)

# If desired, you can write to CSV
# write.csv(kin_df, "lange_kinship_output.csv", row.names=FALSE)
