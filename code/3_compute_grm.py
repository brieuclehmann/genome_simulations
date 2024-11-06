import sys
import tskit
import time
print(tskit.__version__)
print(tskit.__file__)

import pandas as pd
import numpy as np
from scipy.linalg import eigh

sys.path.insert(0, "misc")
import grm_utils
from scipy.sparse.linalg import LinearOperator, eigsh


# Load data
mode = "branch"
path_to_ts = "tree_sequences/chr21/balsac_21_1_sim.ts"
ts = tskit.load(path_to_ts)
selected_df = pd.read_csv("data/balsac_subsample.csv")

# Extract metadata and sample nodes
meta_df = grm_utils.extract_metadata(ts, village_df)
samples_grm = ts.samples()
grm_ind = [i for i in ts.individuals() if int(i.metadata['individual_name']) in list(selected_df['ind'])]

grm_sample_sets = []
for i in grm_ind:
    grm_sample_sets.append([i.nodes[0]])
    grm_sample_sets.append([i.nodes[1]])

start = time.time()
grm = ts.genetic_relatedness_matrix(grm_sample_sets, mode = mode)
end = time.time()
print(end - start)

# Calculate haplotype - individual matrix
n_ind = len(grm_ind)
A = np.zeros((n_ind, 2*n_ind))
for j in range(n_ind):
    A[j, 2*j] = A[j, 2*j + 1] = 1

ind_pca = True
if ind_pca:
    # Calculate haplotype - individual matrix
    n_ind = len(grm_ind)
    A = np.zeros((n_ind, 2*n_ind))
    for j in range(n_ind):
        A[j, 2*j] = A[j, 2*j + 1] = 1
    grm = A @ grm @ A.T

# Save grm as csv
col_names = [str(i.metadata['individual_name']) for i in grm_ind]
grm_df = pd.DataFrame(grm, columns = col_names)
path_to_grm = "output/chr3_grm_" + mode + ".csv"

grm_df.to_csv(path_to_grm, index=False)

n_pc = 10
eigval, eigvec = eigh(grm, subset_by_index=[grm.shape[0] - n_pc, grm.shape[0] - 1])
prin_comp = (grm @ eigvec[:, ::-1]) 


pc_names = ["PC" + str(i+1) for i in range(n_pc)]
pc_df = pd.DataFrame(prin_comp, columns=pc_names)
pc_df['ind'] = col_names
pc_df = pc_df.astype({'ind': 'int64'})

full_df = pd.merge(pc_df, selected_df, on = 'ind')

file_suffix = mode 
full_df.to_csv("output/chr3_pca_" + file_suffix + ".csv", index=False)