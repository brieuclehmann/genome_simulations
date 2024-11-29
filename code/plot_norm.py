import tskit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

tree_type = "recap" # "recap" or ""
ts_recap = tskit.load("tree_sequences/balsac/chr3/sim1/trees_recap.ts")
ts = tskit.load("tree_sequences/balsac/chr3/sim1/trees.ts")

selected_df = pd.read_csv("data/balsac_subsample.csv")
grm_ind = [i for i in ts.individuals() if int(i.metadata['individual_name']) in list(selected_df['ind'])]

# Simplify tree sequence to individuals in grm_ind
sample_nodes = []
for i in grm_ind:
    sample_nodes.append(i.nodes[0])
    sample_nodes.append(i.nodes[1])

ts_sub = ts.simplify(sample_nodes, filter_individuals=False)

# Get sample node IDs for simplified tree sequence
ind_names = [i.metadata['individual_name'] for i in grm_ind]
grm_sub_ind = [i for i in ts_sub.individuals() if i.metadata['individual_name'] in ind_names]
sample_sets = [i.nodes for i in grm_sub_ind]
all_samples = np.array(list({u for s in sample_sets for u in s}))

# Segregating sites / branches
denom_branch = ts.segregating_sites(all_samples, mode = 'branch')
denom_site = ts.segregating_sites(all_samples, mode = 'site')

grm_branch = pd.read_csv("grms/balsac/chr3/sim1/branch_recap.csv")
grm_site = pd.read_csv("grms/balsac/chr3/sim1/site_recap.csv")

grm_branch_norm = grm_branch / denom_branch
grm_site_norm = grm_site / denom_site


# Reorder grms
selected_df = pd.read_csv("data/balsac_subsample.csv")
order_cols = selected_df.sort_values("proband_region")['ind'].tolist()
order_cols = [str(i) for i in order_cols]

def process_grm(grm_df, order_cols):
    grm_df = grm_df.set_index(grm_df.columns.values)
    grm_reorder_df = grm_df[order_cols].reindex(order_cols)
    grm_mat = grm_reorder_df.to_numpy()
    return grm_mat

branch_norm_mat = process_grm(grm_branch_norm, order_cols)
site_norm_mat = process_grm(grm_site_norm, order_cols)

branch_mat = process_grm(grm_branch, order_cols)
site_mat = process_grm(grm_site, order_cols)

# Plot site v branch GRM
plt.clf()
branch_triu = branch_mat[np.triu_indices(branch_mat.shape[0])]
site_triu = site_mat[np.triu_indices(site_mat.shape[0])]
fig, ax = plt.subplots()
ax.scatter(branch_triu, site_triu, s=0.5)
plt.ylabel("Site relatedness")
plt.xlabel("Branch relatedness")
fig.savefig("plots/grm_recap_branch_site_scatter.pdf")


# Plot site v branch GRM
plt.clf()
branch_triu = branch_norm_mat[np.triu_indices(branch_norm_mat.shape[0])]
site_triu = site_norm_mat[np.triu_indices(site_norm_mat.shape[0])]
fig, ax = plt.subplots()
ax.scatter(branch_triu, site_triu, s=0.5)
plt.ylabel("Site relatedness")
plt.xlabel("Branch relatedness")
fig.savefig("plots/grm_recap_branch_site_norm_scatter.pdf")