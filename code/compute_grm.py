import os
import argparse
import sys
import pandas as pd
import numpy as np
import msprime
import time
import demes
import tskit
from scipy.linalg import eigh

sys.path.insert(0, "code")
import pedigree_tools
import grm_tools

def main(args):

    ###################
    #### Load data ####
    ###################

    ts_dir_out = 'tree_sequences/balsac/chr{}/sim{}'.format(
        args.chromosome, args.repetition)

    ts = tskit.load(ts_dir_out + "/trees.ts")
    ts_recap = tskit.load(ts_dir_out + "/trees_recap.ts")

    ######################
    #### Compute GRMs ####
    ######################

    # Prepare output directory
    grm_dir_out = 'grms/balsac/chr{}/sim{}'.format(
        args.chromosome, args.repetition)
    if not os.path.exists(grm_dir_out):
        os.makedirs(grm_dir_out)

    # Load metadata for subsample of 50 individuals from each village
    selected_df = pd.read_csv("data/balsac_subsample.csv")

    # Extract metadata and sample nodes
    grm_ind = [i for i in ts.individuals() if int(i.metadata['individual_name']) in list(selected_df['ind'])]
    col_names = [str(i.metadata['individual_name']) for i in grm_ind]
    for mode in ["site", "branch"]:
        grm = grm_tools.compute_grm(ts, grm_ind, mode = mode)
        grm_df = pd.DataFrame(grm, columns = col_names)
        grm_df.to_csv(grm_dir_out + "/" + mode + ".csv", index=False)
        
        grm_recap = grm_tools.compute_grm(ts_recap, grm_ind, mode = mode)
        grm_recap_df = pd.DataFrame(grm_recap, columns = col_names)
        grm_recap_df.to_csv(grm_dir_out + "/" + mode + "_recap.csv", index=False)

    print("Finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument("-d", "--dir",
        help="directory of input text pedigree"
        )
    parser.add_argument("-o", "--out",
        help="directory of output files"
        )
    parser.add_argument("-ped", "--pedigree_name",
        help="name of input text pedigree (no file extention)"
        )
    # optional arguments
    parser.add_argument("-sfx", "--suffix",
        default="sim",
        help="output file suffix"
        )
    parser.add_argument("-censor", "--censor",
        action="store_true",
        help="removes pedigree information from tree sequence"
        )
    parser.add_argument("-chr", "--chromosome",
        type=int,
        help="specify chromosome number to be simulated"
        )
    parser.add_argument("-m", "--mut_rate",
        default=1.66e-8,
        type=float,
        help="specify mutation rate"
        )
    parser.add_argument("-seed", "--seed_offset",
        default=0,
        type=int,
        help="specify random seed"
        )
    parser.add_argument("-rep", "--repetition",
        default=1,
        type=int,
        help="specify repetition"
        )
    args = parser.parse_args()

    main(args)
