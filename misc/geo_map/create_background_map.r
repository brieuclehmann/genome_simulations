#!/usr/bin/env Rscript

# create_background_map.R
# Builds and saves a background map (bg_map) of Quebec and surrounding areas
# that includes political boundaries, watershed polygons, and simplified watercourses.
# Compatible with modern ggplot2 so as to avoid older internal references (e.g. scales_add_defaults).

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(rmapshaper)
  library(ggplot2)
  library(ggspatial)   # for annotation_scale, annotation_north_arrow
  library(data.table)
})

# ─────────────────────────────────────────────────────────────────────────────
# 1) Define file paths
#    Adjust these to your local environment
# ─────────────────────────────────────────────────────────────────────────────

rdsPath       <- "~/Documents/Genizon/Data/RDS/"    # where your processed RDS files live
#figurePath    <- "~/Documents/Genizon/Genizon_Scripts/Latex/Figures/"
figurePath    <- "~/Documents/genome_simulations_tsrelatedness/misc/geo_map/"
output_map_rds<- file.path(figurePath, "background_map_Jan2025_tsrel.rds")
output_map_jpg<- file.path(figurePath, "background_map_Jan2025_tsrel.jpg")

# Input files
political_file<- file.path(rdsPath, "political_processed.rds")   # simplified political boundaries
wts_file      <- file.path(rdsPath, "wts_processed.rds")         # watershed polygons
water_file    <- file.path(rdsPath, "complete_watercourse_simplified_by_wts.RDS")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Read the processed data
# ─────────────────────────────────────────────────────────────────────────────

message("Reading data...")

political_processed <- readRDS(political_file)
wts_processed       <- readRDS(wts_file)
simplified_water    <- readRDS(water_file)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Identify invalid geometries and remove them if needed
# ─────────────────────────────────────────────────────────────────────────────

message("Checking and fixing invalid geometries...")

# Function to check and report invalid geometries
check_invalid_geometries <- function(sf_obj, name) {
  invalid <- !st_is_valid(sf_obj)
  num_invalid <- sum(invalid, na.rm = TRUE)
  
  if (num_invalid > 0) {
    message(sprintf("WARNING: %d invalid geometries found in %s!", num_invalid, name))
    return(sf_obj[invalid, ])  # Return only the problematic geometries
  } else {
    message(sprintf("✓ All geometries are valid in %s.", name))
    return(NULL)
  }
}

# Run checks and print summary of invalid geometries
invalid_political <- check_invalid_geometries(political_processed, "political_processed")
invalid_wts <- check_invalid_geometries(wts_processed, "wts_processed")
invalid_water <- check_invalid_geometries(simplified_water, "simplified_water")

# Print details of first few invalid geometries (if any)
if (!is.null(invalid_political)) print(head(invalid_political, 3))
if (!is.null(invalid_wts)) print(head(invalid_wts, 3))
if (!is.null(invalid_water)) print(head(invalid_water, 3))

# Fix invalid geometries using `st_make_valid()` where possible
political_processed <- st_make_valid(political_processed)
wts_processed <- st_make_valid(wts_processed)
simplified_water <- st_make_valid(simplified_water)

# Recheck if issues persist
message("Rechecking geometries after st_make_valid()...")

invalid_political <- check_invalid_geometries(political_processed, "political_processed (after fix)")
invalid_wts <- check_invalid_geometries(wts_processed, "wts_processed (after fix)")
invalid_water <- check_invalid_geometries(simplified_water, "simplified_water (after fix)")

# Remove *only* the problematic features if they still exist
if (!is.null(invalid_political)) {
  message("⚠ Removing invalid features from political_processed...")
  political_processed <- political_processed[st_is_valid(political_processed), ]
}

if (!is.null(invalid_wts)) {
  message("⚠ Removing invalid features from wts_processed...")
  wts_processed <- wts_processed[st_is_valid(wts_processed), ]
}

if (!is.null(invalid_water)) {
  message("⚠ Removing invalid features from simplified_water...")
  simplified_water <- simplified_water[st_is_valid(simplified_water), ]
}

# Verify everything is clean
message("Final check completed. Proceeding to the next steps...")

# ─────────────────────────────────────────────────────────────────────────────
# 3) (Optional) Additional filtering or area-based filtering
#    e.g., removing smaller water polygons
# ─────────────────────────────────────────────────────────────────────────────

message("Filtering small polygons or partial features...")

# Example: separate out "Watercourse" from other water to remove small features
# We'll remove water bodies < 1e7 in area
watercourses <- filter(simplified_water, TYPE_TEXT == "Watercourse")
bodies      <- filter(simplified_water, TYPE_TEXT != "Watercourse")

# Compute area, keep only large features
area_info <- data.frame(area = as.numeric(st_area(bodies))) %>%
  mutate(keep = ifelse(area > 1e7, TRUE, FALSE))

bodies_large <- bodies[area_info$keep, ]  # remove small water bodies
water_filtered <- bind_rows(watercourses, bodies_large)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Further simplify or transform the watershed polygons (optional)
# ─────────────────────────────────────────────────────────────────────────────

# e.g., keep them fairly coarse
wts_processed <- wts_processed %>%
  ms_simplify(keep = 0.05)  # e.g. 5% detail

# ─────────────────────────────────────────────────────────────────────────────
# 5) Build the ggplot background map
# ─────────────────────────────────────────────────────────────────────────────

message("Constructing ggplot layers...")

# You can break "political_processed" into subsets, e.g. just Quebec, rest of Canada, etc.
# The code below shows an example, adjusting fill colors accordingly.

bg_map <- ggplot() +
  #  (A) Gray background for Quebec
  geom_sf(
    data = filter(political_processed, juri_en == "Quebec"),
    fill = "gray90", color = NA
  ) +
  #  (B) Slightly different shade for other parts of Canada
  geom_sf(
    data = filter(political_processed, ctry_en == "Canada", juri_en != "Quebec"),
    fill = "gray95", color = NA
  ) +
  #  (C) Even lighter for non-Canadian land
  geom_sf(
    data = filter(political_processed, !ctry_en %in% c("Canada", "Ocean")),
    fill = "gray95", color = NA
  ) +
  #  (D) Watershed polygon outlines in white
  #geom_sf(data = wts_processed, fill = NA, color = "white", size = 0.75) +
  #  (E) Water features as lines or polygons
  geom_sf(data = water_filtered, size = 0.001, color = "gray80", fill = "aliceblue") +
  #  (F) A nice minimal theme with ocean fill
  theme_classic() +
  #  (G) Add scale bar and north arrow
  #annotation_scale(location = "br", width_hint = 0.5) +
  #annotation_north_arrow(
  #  location    = "br", which_north = "true",
  #  pad_x       = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #  style       = north_arrow_fancy_orienteering
  #) +
  theme(
    axis.line        = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "aliceblue")
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6) Save the final map object and optionally a JPG
# ─────────────────────────────────────────────────────────────────────────────

message("Saving bg_map as RDS and JPG...")

saveRDS(bg_map, file = output_map_rds)
ggsave(bg_map, filename = output_map_jpg, width = 24, height = 18, dpi = 300)

message("Done!")
