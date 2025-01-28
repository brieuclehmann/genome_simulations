#!/usr/bin/env Rscript

###############################################################################
# compute_lange_kinship.R
#
# This script computes genealogical (Lange) kinship among probands up to k
# generations using a path-based Breadth-First Search (BFS), then assigns relationship labels.
#
# Usage:
#   Rscript compute_lange_kinship.R pedigree.csv k [output.csv]
#
#   pedigree.csv = CSV with columns: (ind, father, mother)
#                  Father/mother=0 or NA => unknown parent
#   k            = Integer, max number of generations to climb per BFS
#   output.csv   = Optional path to save results. If not provided, results
#                  are printed to console only.
#
# Example:
#   Rscript compute_lange_kinship.R my_pedigree.csv 4 results.csv
#
# Output: A table with columns:
#   proband1, proband2, kinship, relationship
# where:
#   - kinship is the Lange genealogical kinship (sum of (1/2)^(path_len_i+path_len_j))
#     for all shared ancestors within k steps.
#   - relationship is a simplified label (siblings, half_siblings, first_cousins, etc.)
#     based on the minimal path sum for that pair.
#
# Logic Summary:
# 1) For each proband (any individual not appearing as father/mother), run a BFS
#    up to k generations, enumerating every path to an ancestor.
# 2) Self-join these path sets on the common ancestor column to pair up
#    (proband1->ancestor path) with (proband2->ancestor path).
# 3) Each pair of paths contributes (1/2)^(path_len1 + path_len2) to kinship.
# 4) Sum those contributions over all shared ancestors => Lange kinship.
# 5) For labeling, look at the minimal path sum among all matched paths, then
#    define relationships (2 => siblings, 4 => first_cousins, 6 => second_cousins, etc.).
# 6) For siblings, check father/mother IDs to see if they're "full" or "half".
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

###############################################################################
# 1) Enumerate all ancestor paths (BFS, up to k)
###############################################################################

# This function performs a BFS from `focal_ind` upward (father, mother)
# collecting all "paths" up to k generations. Each path is stored as a row in
# the returned data frame. The BFS approach ensures we capture multiple routes
# to the same ancestor if inbreeding or complex loops exist.
#
# Columns in the returned data:
#   proband     = the focal individual ID
#   ancestor    = an ancestor ID reached by BFS
#   path_length = number of edges (parent->child steps) from proband to ancestor
#   path_string = textual representation of the path 

enumerate_paths_up_to_k <- function(pedigree, focal_ind, k) {
  father_map <- setNames(pedigree$father, pedigree$ind)
  mother_map <- setNames(pedigree$mother, pedigree$ind)

  results <- data.frame(
    proband     = integer(),
    ancestor    = integer(),
    path_length = integer(),
    path_string = character(),
    stringsAsFactors = FALSE
  )

  # Our queue is a list of BFS "states". Each state has:
  #   - current_id : current ancestor ID
  #   - length     : how many edges from the focal_ind
  #   - path_str   : textual chain of IDs
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

    # Record the current BFS state as a row in 'results'
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

# A helper that enumerates paths for multiple probands and stacks them together.
collect_paths_for_probands <- function(pedigree, probands, k) {
  do.call(rbind, lapply(probands, function(p) {
    enumerate_paths_up_to_k(pedigree, p, k)
  }))
}

###############################################################################
# 2) Compute Lange kinship by summing (1/2)^(path_len1 + path_len2)
###############################################################################

# This function:
# 1) Gathers all BFS paths from each proband
# 2) Self-joins them on 'ancestor' to find pairs of proband paths that converge
#    at the same ancestor
# 3) Each pair of paths contributes (1/2)^(path_len1 + path_len2)
# 4) Sums the result to yield total genealogical kinship
#
# Return: list with
#   - kinship: a data frame of (pair_id, proband1, proband2, kinship)
#   - joined_paths: the full row-level join (for second-pass labeling)

compute_lange_kinship_up_to_k <- function(pedigree, probands, k) {
  # 1) BFS for all probands
  all_paths <- collect_paths_for_probands(pedigree, probands, k)

  # 2) We'll label them as paths_l and paths_r, then do an "inner_join" on ancestor
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

  # relationship="many-to-many" means we allow multiple BFS rows to match
  # multiple BFS rows (e.g. in inbreeding or multi-route scenarios).
  joined <- paths_l %>%
    inner_join(paths_r, by = "ancestor", relationship = "many-to-many") %>%
    # Ensure we don't double-count symmetrical pairs (p1,p2) vs (p2,p1)
    filter(proband1 < proband2) %>%
    mutate(
      pair_id = paste0(proband1, "_", proband2),
      contrib = (1/2)^(path_len1 + path_len2)
    )

  # 3) Sum up contributions by pair
  kin_df <- joined %>%
    group_by(pair_id, proband1, proband2) %>%
    summarise(
      kinship = sum(contrib),
      .groups = "drop"
    ) %>%
    arrange(proband1, proband2)

  # Return both the aggregated kinship and the row-level joined table
  list(kinship = kin_df, joined_paths = joined)
}

###############################################################################
# 3) Relationship labeling based on minimal path sum
###############################################################################
#
# We'll find, for each pair (p1,p2), the minimal total path length = (path_len1 + path_len2)
# across all their shared-ancestor paths. That minimal sum usually defines
# the "closest relationship" genealogically. For example:
#   - sum=2 => siblings (one parent each)
#   - sum=4 => first cousins
#   - sum=6 => second cousins
#   - sum=2n => (n-1)_degree_relatives for higher n
# Then we do a specialized check for siblings to see if they share the
# same father and mother => full_siblings or half_siblings.

label_relationships <- function(joined_paths, pedigree) {
  joined_paths <- joined_paths %>%
    mutate(total_path_len = path_len1 + path_len2)

  # 1) Find minimal sum per pair
  min_path_summary <- joined_paths %>%
    group_by(proband1, proband2) %>%
    summarise(
      min_sum = min(total_path_len),
      .groups = "drop"
    )

  # 2) Define a function that maps min_sum to a base relationship label
  #    We skip half-cousins. For siblings, we do further checks.
  classify_base_relationship <- function(s) {
    if (s == 2) {
      return("siblings")
    } else if (s == 4) {
      return("first_cousins")
    } else if (s == 6) {
      return("second_cousins")
    } else if (s > 2 && (s %% 2) == 0) {
      # For sum >=8, e.g. 8 => 3_degree_relatives, etc.
      g <- (s / 2) - 1
      return(paste0(g, "_degree_relatives"))
    } else {
      return("unclassified")
    }
  }

  # 3) For siblings, check father/mother for full vs half
  label_sibling_type <- function(p1, p2, ped) {
    f1 <- ped$father[ped$ind == p1]
    m1 <- ped$mother[ped$ind == p1]
    f2 <- ped$father[ped$ind == p2]
    m2 <- ped$mother[ped$ind == p2]

    father_match <- (!is.na(f1) && !is.na(f2) && f1 == f2 && f1 != 0)
    mother_match <- (!is.na(m1) && !is.na(m2) && m1 == m2 && m1 != 0)

    if (father_match && mother_match) {
      return("full_siblings")
    } else if (father_match || mother_match) {
      return("half_siblings")
    } else {
      # if father/mother don't match at all, we call them "siblings" anyway
      # (maybe data is incomplete, or adoptive scenario).
      return("siblings")
    }
  }

  # 4) Assign the base relationship. Then refine siblings => full/half_siblings
  min_path_summary <- min_path_summary %>%
    mutate(
      base_relationship = sapply(min_sum, classify_base_relationship),
      final_relationship = base_relationship
    ) %>%
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

# Command line interface:
#   Rscript compute_lange_kinship.R <pedigree.csv> <k> [<output.csv>]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript compute_lange_kinship.R pedigree.csv k [output.csv]\n")
  quit(status = 1)
}

ped_file <- args[1]
k <- as.integer(args[2])
output_csv <- if (length(args) >= 3) args[3] else NULL

# Read pedigree from CSV
pedigree <- read.csv(ped_file, stringsAsFactors = FALSE)

# Convert father=0 or mother=0 => NA for unknown
pedigree <- pedigree %>%
  mutate(
    father = ifelse(father == 0, NA, father),
    mother = ifelse(mother == 0, NA, mother)
  )

# Define "probands" as those not appearing as father or mother
all_parents <- c(pedigree$father, pedigree$mother)
probands <- setdiff(pedigree$ind, all_parents[!is.na(all_parents)])

# 1) Compute Lange kinship from enumerated BFS paths
res <- compute_lange_kinship_up_to_k(pedigree, probands, k)
kin_df <- res$kinship
joined_paths <- res$joined_paths

# 2) Label relationships (siblings, cousins, etc.) based on minimal path sums
rel_df <- label_relationships(joined_paths, pedigree)

# 3) Merge kinship with final relationship labels
out_df <- kin_df %>%
  left_join(
    rel_df %>% select(proband1, proband2, final_relationship),
    by = c("proband1", "proband2")
  ) %>%
  select(proband1, proband2, kinship, relationship = final_relationship) %>%
  arrange(proband1, proband2)

cat(sprintf("Lange kinship among probands (up to %d generations) + relationship labels:\n", k))
print(out_df %>% arrange(kinship), n = 10)

# Optionally write to output CSV
if (!is.null(output_csv)) {
  write.csv(out_df, output_csv, row.names = FALSE)
  cat(sprintf("Wrote results to %s\n", output_csv))
}
