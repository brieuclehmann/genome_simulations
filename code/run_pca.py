import os
import argparse
import sys
import time
sys.path.insert(0,"../tskit/python")
import tskit
print(tskit.__version__)

import pandas as pd
import numpy as np
import stdpopsim    
print(stdpopsim.__version__)

from scipy.linalg import eigh
import matplotlib
import matplotlib.pyplot as plt

sys.path.insert(0, "code")
import pedigree_tools
import grm_tools

def main(args):

    ###################
    #### Load data ####
    ###################

    ts_dir_out = 'tree_sequences/balsac/chr3/sim1'

    ts = tskit.load(ts_dir_out + "/trees_recap.ts")

    ### Rescale tree sequence ###
    ts_path = ts_dir_out + "/trees_recap_rescaled.ts"
    if os.path.exists(ts_path):
        new_ts = tskit.load(ts_path)
    else:
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr3", genetic_map="HapMapII_GRCh37")
        rate_map = contig.recombination_map
        t = ts.dump_tables()
        t.mutations.clear()
        t.sites.clear()
        new_left = rate_map.get_cumulative_mass(t.edges.left)
        new_right = rate_map.get_cumulative_mass(t.edges.right)
        t.edges.set_columns(left=new_left, right=new_right, parent=t.edges.parent, child=t.edges.child)
        new_ts = t.tree_sequence()
        new_ts.dump(ts_path)

    ######################
    #### Load GRMs ####
    ######################

    # Prepare output directory
    grm_dir_out = 'grms/balsac/chr3/sim1'

    relatives_df = pd.read_csv("data/balsac_relatives.csv")
    pedigree_df = pd.read_csv("data/ascendance.csv", sep=";")
    region_df = pedigree_df[['FiveRegions','ProbandID']].drop_duplicates()
    region_df = region_df.rename(columns={'FiveRegions':'proband_region','ProbandID':'ind'})
    full_df = pd.read_csv("data/balsac_pedigree.csv")
    proband_df = full_df[~(full_df['ind'].isin(full_df['mother'])) & ~(full_df['ind'].isin(full_df['father']))]
    proband_df = proband_df.merge(region_df, on='ind')

    n_select = 20
    selected_df = proband_df.groupby('proband_region').sample(n_select, random_state=42)
    selected_df = proband_df
    grm_ind = [i for i in ts.individuals() if int(i.metadata['individual_name']) in list(selected_df['ind'])]
    col_names = [str(i.metadata['individual_name']) for i in grm_ind]

    # Compute GRM plus rescaled GRM
    mode = 'branch'
    file_path = grm_dir_out + "/" + mode + "pca.csv"
    if os.path.exists(file_path):
        grm_df = pd.read_csv(file_path)
        grm = grm_df.to_numpy()
    else:
        grm = grm_tools.compute_grm(ts, grm_ind, mode = mode)
        grm_df = pd.DataFrame(grm, columns = col_names)
        grm_df.to_csv(grm_dir_out + "/" + mode + "pca.csv", index=False)

    # Recaled GRM
    file_path = grm_dir_out + "/" + mode + "pca_rescaled.csv"
    if os.path.exists(file_path):
        grm_rescaled_df = pd.read_csv(file_path)
        grm_rescaled = grm_rescaled_df.to_numpy()
    else:
        grm_rescaled = grm_tools.compute_grm(new_ts, grm_ind, mode = mode)
        grm_rescaled_df = pd.DataFrame(grm_rescaled, columns = col_names)
        grm_rescaled_df.to_csv(grm_dir_out + "/" + mode + "pca_rescaled.csv", index=False)


    # Extract metadata and sample nodes
    selected_ind = [i for i in new_ts.individuals() if int(i.metadata['individual_name']) in list(selected_df['ind'])]
    pca_ind = [i.id for i in selected_ind]
    balsac_ind = [int(i.metadata['individual_name']) for i in selected_ind]

    label_df = pd.DataFrame({"ind":balsac_ind})
    # label_df = pd.DataFrame({"ind":grm.columns.values})
    label_df = label_df.astype({"ind":np.int32})
    label_df = label_df.merge(region_df, on = 'ind', how = 'left')
    labels = label_df['proband_region']
    labels_int, labels_unique = pd.factorize(labels)

    col_names = [str(i.metadata['individual_name']) for i in grm_ind]

    #####################
    #### Compute PCs ####
    #####################

    n_pc = 10
    pc_names = ["PC" + str(i+1) for i in range(n_pc)]

    def plot_pca(prin_comp, file_path):
        colors = matplotlib.colormaps.get_cmap('Set2')(np.arange(5))
        fig, ax = plt.subplots(figsize=(5.2,5))
        for i, label in enumerate(labels_unique):
            idx = labels_int == i
            ax.scatter(prin_comp[idx,0], prin_comp[idx,1], s=5, label=label)

        ax.legend(frameon=False)
        plt.show()
        plt.savefig(file_path)
    

    def plot_pca(prin_comp, labels_int, labels_unique, file_path):
        """
        Plots the first 4 pairs of PC plots (PC1 vs PC2, PC3 vs PC4, etc.).

        Parameters:
            prin_comp (ndarray): Principal components array of shape (n_samples, n_components).
            labels_int (ndarray): Integer labels for data points.
            labels_unique (list): Unique labels corresponding to the integer labels.
            file_path (str): Path to save the final plot.
        """
        colors = plt.cm.Set2(np.arange(len(labels_unique)))  # Define a colormap for the labels
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))    # Create 2x2 subplots
        axes = axes.flatten()  # Flatten the axes array for easier indexing

        pc_pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]  # Define the first 4 PC pairs

        for idx, (pc_x, pc_y) in enumerate(pc_pairs):
            ax = axes[idx]
            for i, label in enumerate(labels_unique):
                mask = labels_int == i  # Filter points by label
                ax.scatter(prin_comp[mask, pc_x], prin_comp[mask, pc_y], s=5, label=label, color=colors[i])
            
            ax.set_title(f'PC{pc_x + 1} vs PC{pc_y + 1}')
            ax.set_xlabel(f'PC{pc_x + 1}')
            ax.set_ylabel(f'PC{pc_y + 1}')
            ax.legend(frameon=False, fontsize=8)

        plt.tight_layout()
        plt.savefig(file_path)  # Save the figure

# Example call
# plot_pca(prin_comp, labels_int, labels_unique, 'pca_plots.png')



    # eigh PCA
    # eigval, eigvec = eigh(grm, subset_by_index=[grm.shape[0] - n_pc, grm.shape[0] - 1])
    # prin_comp = (grm @ eigvec[:, ::-1]) 
    # file_path = "plots/test_pca_grm.png"
    # plot_pca(prin_comp, file_path)

    # # eigh PCA rescaled
    # eigval, eigvec = eigh(grm_rescaled, subset_by_index=[grm_rescaled.shape[0] - n_pc, grm_rescaled.shape[0] - 1])
    # prin_comp = (grm_rescaled @ eigvec[:, ::-1]) 
    # file_path = "plots/test_pca_grm_rescaled.png"
    # plot_pca(prin_comp, file_path)

    # ts PCA
    out = ts.pca(individuals=np.asarray(pca_ind), iterated_power=5, random_seed=1, num_components=10)
    prin_comp = out.factors
    file_path = "plots/test_pca_ts.png"
    plot_pca(prin_comp, labels_int, labels_unique, file_path)

    # ts rescaled PCA
    start_time = time.time()
    out = new_ts.pca(individuals=np.asarray(pca_ind), iterated_power=5, random_seed=1, num_components=10)
    end_time = time.time()
    print(end_time - start_time)
    prin_comp = out.factors
    file_path = "plots/test_pca_ts_rescaled.png"
    plot_pca(prin_comp, labels_int, labels_unique, file_path)


    # pc_df = pd.DataFrame(prin_comp, columns=pc_names)
    # pc_df['ind'] = col_names
    # pc_df = pc_df.astype({'ind': 'int64'})

    # # Save PCs
    # pca_dir_out = 'pca/balsac/chr{}/sim{}'.format(
    #        args.chromosome, args.repetition)
    # if not os.path.exists(pca_dir_out):
    #     os.makedirs(pca_dir_out)
    # full_df = pd.merge(pc_df, selected_df, on = 'ind')
    # full_df.to_csv(pca_dir_out + '/pca.csv', index=False)

    # # Recapitated
    # eigval, eigvec = eigh(grm_recap, subset_by_index=[grm_recap.shape[0] - n_pc, grm_recap.shape[0] - 1])
    # prin_comp = (grm_recap @ eigvec[:, ::-1]) 

    # pc_df = pd.DataFrame(prin_comp, columns=pc_names)
    # pc_df['ind'] = col_names
    # pc_df = pc_df.astype({'ind': 'int64'})

    # full_df = pd.merge(pc_df, selected_df, on = 'ind')
    # full_df.to_csv(pca_dir_out + '/pca_recap.csv', index=False)

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
