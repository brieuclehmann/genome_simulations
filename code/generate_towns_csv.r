#!/usr/bin/env Rscript

# generate_towns_csv.R
# 1) Defines a table of parishes of interest (region, name, lieu).
# 2) Reads the BALSAC metadata (which includes 'lieu', 'Lat', 'Lon', etc.).
# 3) Joins them, then writes 'ts-relatedness-towns.csv'.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# Hard-coded paths (edit as needed)
meta_file <- "/Users/luke/Documents/Genizon/BALSAC/Balsac_aout_2021_v2/balsac_metadata.csv"
output_csv <- "/Users/luke/Desktop/ts-relatedness-towns.csv"

# 1) Define the parishes of interest
towns_of_interest <- data.frame(
  region = c(
    "Chaleur Bay","Chaleur Bay","Chaleur Bay","Chaleur Bay","Chaleur Bay",
    "Batiscan","Batiscan","Batiscan","Batiscan",
    "Chaudière","Chaudière","Chaudière","Chaudière","Chaudière",
    "L’Assomption","L’Assomption","L’Assomption","L’Assomption",
    "Mistassini","Mistassini","Mistassini","Mistassini","Mistassini"
  ),
  name = c(
    "St Michel","St François De Sales","St Georges De Malbaie","St Pierre De Malbaie","St Joseph",
    "Ste Geneviève De Batiscan","St Luc De Vincennes","St Narcisse","St Stanislas",
    "St Georges","St Benoit Labre","St Philibert","St Come","St Martin De Tours",
    "St Jacques","St Alexis","Ste Marie Salomée","St Esprit",
    "St Félicien","St Méthode","Notre Dame De La Dore","St Cyrille","Ste Lucie"
  ),
  lieu = c(
    11540,12187,11223,12432,13601,
    11143,11045,11057,11648,
    13415,10897,10913,10981,11497,
    11017,11349,11151,10991,
    10731,2984,11235,11578,11333
  ),
  stringsAsFactors = FALSE
)

# 2) Read the BALSAC metadata
meta <- fread(meta_file)  # Must contain at least: lieu, Lat, Lon, etc.

# 3) Merge the two sets on "lieu" and select relevant columns
parish_data <- towns_of_interest %>%
  left_join(meta, by = "lieu") %>%
  # Keep whichever columns you want in the final CSV
  select(region, name, t, lieu, Lat, Lon)

# 4) Write to CSV
fwrite(parish_data, file = output_csv)
cat("Wrote CSV to:", output_csv, "\n")
