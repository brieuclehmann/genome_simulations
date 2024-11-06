import sys
import tskit
import time
print(tskit.__version__)
print(tskit.__file__)

import pandas as pd
import numpy as np

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

# Calculate haplotype - individual matrix
n_ind = len(grm_ind)
A = np.zeros((n_ind, 2*n_ind))
for j in range(n_ind):
    A[j, 2*j] = A[j, 2*j + 1] = 1

ind_pca = True
if ind_pca:
    grm_matvec = lambda x: A @ ts.genetic_relatedness_vector(x @ A, mode = "branch")
    mv_linop = LinearOperator((n_ind, n_ind), matvec=grm_matvec)
else:
    grm_matvec = lambda x: ts.genetic_relatedness_vector(x, mode = "branch")
    mv_linop = LinearOperator((2 * n_ind, 2 * n_ind), matvec=grm_matvec)

n_pc = 10
start = time.time()
eigval, eigvec = eigsh(mv_linop, n_pc)
end = time.time()
print(end - start)


prin_comp = (mv_linop @ eigvec[:, ::-1]) 

row_names = [i.id for i in grm_sample_ind]
col_names = ["PC" + str(i+1) for i in range(n_pc)]
pc_df = pd.DataFrame(prin_comp, columns=col_names)

grm_id = [int(i) for i in grm_df.columns.values]
village_names = meta_df.loc[meta_df['id'].isin(row_names)]['Nom']
pc_df['Village'] = list(village_names)

file_suffix = mode + "_matvec_prop"
if snp_keep < 1:
    file_suffix = mode + "_matvec_prop" + str(int(snp_keep * 100))

pc_df.to_csv("output/chr3_pca_" + file_suffix + ".csv", index=False)