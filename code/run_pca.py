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
import matplotlib.pyplot as plt

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
    #### Load GRMs ####
    ######################

    # Prepare output directory
    grm_dir_out = 'grms/balsac/chr{}/sim{}'.format(
        args.chromosome, args.repetition)
    if not os.path.exists(grm_dir_out):
        os.makedirs(grm_dir_out)

    mode = "branch"
    grm = pd.read_csv(grm_dir_out + "/" + mode + ".csv")

    # Load metadata for subsample of 100 individuals from each village
    selected_df = pd.read_csv("data/balsac_subsample.csv")
    region_df  = selected_df[['ind' , 'proband_region']]
    region_df['ind'] = region_df['ind'].astype(np.int32)

    label_df = pd.DataFrame({"ind":grm.columns.values})
    label_df['ind'] = label_df['ind'].astype(np.int32)
    label_df = label_df.merge(region_df, on = 'ind', how = 'left')
    labels = label_df['proband_region']
    labels_int, labels_unique = pd.factorize(labels)

    # Extract metadata and sample nodes
    grm_ind = [i for i in ts.individuals() if int(i.metadata['individual_name']) in list(selected_df['ind'])]
    col_names = [str(i.metadata['individual_name']) for i in grm_ind]

    #####################
    #### Compute PCs ####
    #####################

    n_pc = 10
    pc_names = ["PC" + str(i+1) for i in range(n_pc)]

    grm = grm.to_numpy()
    eigval, eigvec = eigh(grm, subset_by_index=[grm.shape[0] - n_pc, grm.shape[0] - 1])
    prin_comp = (grm @ eigvec[:, ::-1]) 

    colors = matplotlib.colormaps.get_cmap('Set2')(np.arange(5))
    fig, ax = plt.subplots(figsize=(5.2,5))
    for i, label in enumerate(labels_unique):
        idx = labels_int == i
        ax.scatter(prin_comp[idx,2], prin_comp[idx,3], s=5, label=label)

    ax.legend(frameon=False)
    plt.show()
    plt.savefig("plots/test.png")

    pc_df = pd.DataFrame(prin_comp, columns=pc_names)
    pc_df['ind'] = col_names
    pc_df = pc_df.astype({'ind': 'int64'})

    # Save PCs
    pca_dir_out = 'pca/balsac/chr{}/sim{}'.format(
           args.chromosome, args.repetition)
    if not os.path.exists(pca_dir_out):
        os.makedirs(pca_dir_out)
    full_df = pd.merge(pc_df, selected_df, on = 'ind')
    full_df.to_csv(pca_dir_out + '/pca.csv', index=False)

    # Recapitated
    eigval, eigvec = eigh(grm_recap, subset_by_index=[grm_recap.shape[0] - n_pc, grm_recap.shape[0] - 1])
    prin_comp = (grm_recap @ eigvec[:, ::-1]) 

    pc_df = pd.DataFrame(prin_comp, columns=pc_names)
    pc_df['ind'] = col_names
    pc_df = pc_df.astype({'ind': 'int64'})

    full_df = pd.merge(pc_df, selected_df, on = 'ind')
    full_df.to_csv(pca_dir_out + '/pca_recap.csv', index=False)

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
