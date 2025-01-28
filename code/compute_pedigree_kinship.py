#!/usr/bin/env python3

import sys
import pandas as pd
import xarray as xr
import numpy as np

# sgkit is needed:
#    pip install sgkit
import sgkit


def pedigree_to_sgkit_dataset(pedigree_df):
    """
    Convert a pedigree DataFrame with columns ['ind', 'father', 'mother']
    to an xarray Dataset suitable for sgkit.pedigree_kinship.
    - Expects integer IDs for ind, father, mother.
    - father=0 or mother=0 indicates no parent.
    """

    # Sort or otherwise ensure the data frame is in a stable order
    pedigree_df = pedigree_df.sort_values("ind").reset_index(drop=True)

    # Create sample_id strings (e.g. "S1", "S2", ...)
    sample_ids = [f"S{int(i)}" for i in pedigree_df["ind"]]

    # Build a map: integer ID -> "S#" string
    # If father or mother is 0, map that to "."
    id_map = {row["ind"]: f"S{int(row['ind'])}" for _, row in pedigree_df.iterrows()}
    # For missing parents
    id_map[0] = "."

    # Create parent_id array (2 columns => father, mother)
    # shape: (n_samples, n_parents)
    n_samples = pedigree_df.shape[0]
    parent_id = np.full((n_samples, 2), ".", dtype=object)

    # father is in column 0, mother is in column 1
    father_ids = [id_map.get(fa, ".") for fa in pedigree_df["father"]]
    mother_ids = [id_map.get(mo, ".") for mo in pedigree_df["mother"]]
    parent_id[:, 0] = father_ids
    parent_id[:, 1] = mother_ids

    # Build xarray.Dataset
    ds = xr.Dataset()
    ds = ds.assign_coords(samples=("samples", np.arange(n_samples)))
    ds["sample_id"] = ("samples", sample_ids)

    # We use dimension names ("samples", "parents") for parent_id
    ds["parent_id"] = (("samples", "parents"), parent_id)

    return ds


def compute_pedigree_kinship(input_csv, output_csv=None):
    """
    Reads the pedigree CSV (with at least 'ind', 'father', 'mother' columns),
    computes pairwise kinship using sgkit, and writes a long-form CSV of
    pairwise kinship with columns: sample_i, sample_j, kinship.
    """

    # 1. Read the pedigree from CSV
    df = pd.read_csv(input_csv)
    # Ensure father/mother are numeric
    df["father"] = pd.to_numeric(df["father"], errors="coerce").fillna(0).astype(int)
    df["mother"] = pd.to_numeric(df["mother"], errors="coerce").fillna(0).astype(int)

    # 2. Convert to an sgkit Dataset
    ds = pedigree_to_sgkit_dataset(df)

    # 3. Run pedigree_kinship
    #    - This uses the default 'diploid' method
    #    - If you have autopolyploids, you would add `method="Hamilton-Kerr"`
    ds_kinship = sgkit.pedigree_kinship(ds, method="diploid", merge=True)

    # The result has ds_kinship["stat_pedigree_kinship"] as a 2D matrix
    kinship_matrix = ds_kinship["stat_pedigree_kinship"].values
    sample_ids = ds_kinship["sample_id"].values

    # 4. Convert the kinship matrix to a long table
    #    with columns sample_i, sample_j, kinship
    n = len(sample_ids)
    rows = []
    for i in range(n):
        for j in range(n):
            rows.append((sample_ids[i], sample_ids[j], kinship_matrix[i, j]))

    kinship_df = pd.DataFrame(rows, columns=["sample_i", "sample_j", "kinship"])

    # Optionally, write to CSV or return the DataFrame
    if output_csv:
        kinship_df.to_csv(output_csv, index=False)
        print(f"Kinship matrix saved to {output_csv}")
    else:
        print(kinship_df.head(20))

    return kinship_df


if __name__ == "__main__":
    """
    Example command-line usage:
        python compute_pedigree_kinship.py input_pedigree.csv output_kinship.csv

    input_pedigree.csv must have columns:
        ind,father,mother, (optionally sex, or other columns are ignored)
    """
    if len(sys.argv) < 2:
        print("Usage: python compute_pedigree_kinship.py input_pedigree.csv [output_kinship.csv]")
        sys.exit(1)

    input_csv = sys.argv[1]
    if len(sys.argv) > 2:
        output_csv = sys.argv[2]
    else:
        output_csv = None

    compute_pedigree_kinship(input_csv, output_csv)
