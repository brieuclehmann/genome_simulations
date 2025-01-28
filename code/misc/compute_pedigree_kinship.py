#!/usr/bin/env python3

import sys
import pandas as pd
import xarray as xr
import numpy as np

# sgkit is needed:
#    pip install sgkit
import sgkit


def truncate_pedigree_at_depth(df, k, probands=None):
    """
    Truncate the pedigree so that ancestors more than k generations above
    the probands are treated as founders (father=0, mother=0).

    Parameters
    ----------
    df : pd.DataFrame
        Pedigree table with columns 'ind', 'father', 'mother' (all integers).
        0 indicates unknown parent.
    k : int
        Maximum number of generations to keep above each proband.
    probands : list[int] or None
        List of individual IDs to start from. If None, probands are
        automatically defined as individuals not appearing in any father or
        mother column.

    Returns
    -------
    pd.DataFrame
        A copy of the pedigree dataframe, where individuals deeper than k
        generations from the probands have father=0 and mother=0.
    """

    df = df.copy()

    # 1) If probands are not specified, define them:
    #    probands = all individuals who are not listed as a father/mother
    if probands is None:
        all_fathers = set(df["father"]) - {0}
        all_mothers = set(df["mother"]) - {0}
        parent_set = all_fathers.union(all_mothers)
        all_inds = set(df["ind"])
        probands = list(all_inds - parent_set)

    # 2) Build quick father/mother lookups
    father_map = {}
    mother_map = {}
    for row in df.itertuples(index=False):
        # row is: Pandas(ind=..., father=..., mother=...)
        father_map[row.ind] = row.father
        mother_map[row.ind] = row.mother

    # 3) BFS/DFS from each proband, tracking depth
    from collections import deque

    depth_map = {}  # depth_map[ind] = generations above a proband
    queue = deque()

    # Initialize the queue with each proband at depth 0
    for p in probands:
        depth_map[p] = 0
        queue.append(p)

    while queue:
        child = queue.popleft()
        d = depth_map[child]

        if d >= k:
            # No need to go further up from here
            continue

        # Identify the father and mother
        fa = father_map.get(child, 0)
        mo = mother_map.get(child, 0)

        # For each nonzero parent, if not visited, set depth = d+1
        for parent_id in (fa, mo):
            if parent_id != 0 and parent_id not in depth_map:
                depth_map[parent_id] = d + 1
                queue.append(parent_id)

    # 4) For anyone who is deeper than k generations (i.e., depth > k)
    #    or not in depth_map at all, set father/mother = 0 (founder).
    #    Also if the father or mother is not in depth_map or depth_map > k.
    for idx, row in df.iterrows():
        fa, mo = row["father"], row["mother"]
        # If father is outside the k-depth
        if fa != 0:
            if fa not in depth_map or depth_map[fa] > k:
                df.at[idx, "father"] = 0
        # If mother is outside the k-depth
        if mo != 0:
            if mo not in depth_map or depth_map[mo] > k:
                df.at[idx, "mother"] = 0

    return df


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
    n_samples = pedigree_df.shape[0]
    parent_id = np.full((n_samples, 2), ".", dtype=object)

    father_ids = [id_map.get(fa, ".") for fa in pedigree_df["father"]]
    mother_ids = [id_map.get(mo, ".") for mo in pedigree_df["mother"]]
    parent_id[:, 0] = father_ids
    parent_id[:, 1] = mother_ids

    # Build xarray.Dataset
    ds = xr.Dataset()
    ds = ds.assign_coords(samples=("samples", np.arange(n_samples)))
    ds["sample_id"] = ("samples", sample_ids)

    ds["parent_id"] = (("samples", "parents"), parent_id)

    return ds


def compute_pedigree_kinship(input_csv, output_csv=None, max_depth=None):
    """
    Reads the pedigree CSV (with at least 'ind', 'father', 'mother' columns),
    optionally truncates it to max_depth generations, computes pairwise kinship
    using sgkit, and writes a long-form CSV of pairwise kinship with columns:
    sample_i, sample_j, kinship.

    Parameters
    ----------
    input_csv : str
        Path to input pedigree CSV.
    output_csv : str or None
        If provided, writes the resulting kinship table to CSV.
    max_depth : int or None
        If an integer, then pedigree is truncated to max_depth generations
        above each proband (individual not listed as mother/father).
    """

    # 1. Read the pedigree from CSV
    df = pd.read_csv(input_csv)
    # Ensure father/mother are numeric
    df["father"] = pd.to_numeric(df["father"], errors="coerce").fillna(0).astype(int)
    df["mother"] = pd.to_numeric(df["mother"], errors="coerce").fillna(0).astype(int)

    # 2. Truncate if requested
    if max_depth is not None:
        df = truncate_pedigree_at_depth(df, k=max_depth)

    # 3. Convert to an sgkit Dataset
    ds = pedigree_to_sgkit_dataset(df)

    # 4. Run pedigree_kinship
    ds_kinship = sgkit.pedigree_kinship(
        ds,
        method="diploid",
        allow_half_founders=True,
        merge=True
    )

    kinship_matrix = ds_kinship["stat_pedigree_kinship"].values
    sample_ids = ds_kinship["sample_id"].values

    # 5. Convert the kinship matrix to a long table: (sample_i, sample_j, kinship)
    n = len(sample_ids)
    rows = []
    for i in range(n):
        for j in range(n):
            rows.append((sample_ids[i], sample_ids[j], kinship_matrix[i, j]))

    kinship_df = pd.DataFrame(rows, columns=["sample_i", "sample_j", "kinship"])

    # 6. Output or print the kinship table
    if output_csv:
        kinship_df.to_csv(output_csv, index=False)
        print(f"Kinship matrix saved to {output_csv}")
    else:
        print(kinship_df.head(20))

    return kinship_df


if __name__ == "__main__":
    """
    Example usage at the command line:

        python compute_pedigree_kinship.py input_pedigree.csv output_kinship.csv 5

    Here, 5 is the maximum depth to climb above each proband. If you omit it,
    the script calculates kinship with no depth limit.
    """
    if len(sys.argv) < 2:
        print(
            "Usage: python compute_pedigree_kinship.py input_pedigree.csv [output_kinship.csv] [max_depth]"
        )
        sys.exit(1)

    input_csv = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else None
    max_depth = int(sys.argv[3]) if len(sys.argv) > 3 else None

    compute_pedigree_kinship(input_csv, output_csv, max_depth)
