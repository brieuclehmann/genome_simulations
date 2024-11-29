import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Get selected individuals in region order
selected_df = pd.read_csv("data/balsac_subsample.csv")
order_cols = selected_df.sort_values("proband_region")['ind'].tolist()
order_cols = [str(i) for i in order_cols]

def process_grm(grm_df, order_cols):
    grm_df = grm_df.set_index(grm_df.columns.values)
    grm_reorder_df = grm_df[order_cols].reindex(order_cols)
    grm_mat = grm_reorder_df.to_numpy()
    return grm_mat

# Load in data
grm_file = "grms/balsac/chr3/sim1/branch.csv"
grm_recap_file = "grms/balsac/chr3/sim1/branch_recap.csv"
prm_file = "output/chr3_prm.csv"

grm_df = pd.read_csv(grm_file)
grm_recap_df = pd.read_csv(grm_recap_file)

prm_df = pd.read_csv(prm_file)
ped_df = pd.read_csv("data/balsac_pedigree.csv")
ind_dict = dict(zip(ped_df.pedID, ped_df.ind))
new_cols = [str(ind_dict[int(i)]) for i in prm_df.columns.values]
prm_df.columns = new_cols

# Process GRMs (i.e. reorder individuals)
prm_mat = process_grm(prm_df, order_cols)
grm_mat = process_grm(grm_df, order_cols)
grm_recap_mat = process_grm(grm_recap_df, order_cols)

relative_df = pd.read_csv("data/balsac_relatives.csv")
relative_df = relative_df.loc[
    (relative_df['proband1'].astype("string").isin(new_cols)) &
    (relative_df['proband2'].astype("string").isin(new_cols)) 
    ]

region_df  = selected_df[['ind' , 'proband_region']]
region_df['ind'] = region_df['ind'].astype(np.int32)


# Genetic relatedness v pedigree relatedness
plt.clf()
grm_triu = grm_mat[np.triu_indices(grm_mat.shape[0])]
prm_triu = prm_mat[np.triu_indices(prm_mat.shape[0])]
fig, ax = plt.subplots()
ax.scatter(prm_triu, grm_triu)
fig.savefig("plots/grm_prm_scatter.pdf")

# GRM v recap-GRM
plt.clf()
grm_triu = grm_mat[np.triu_indices(grm_mat.shape[0])]
grm_recap_triu = grm_recap_mat[np.triu_indices(grm_recap_mat.shape[0])]
fig, ax = plt.subplots()
ax.scatter(grm_recap_triu, grm_triu, s=0.5)
fig.savefig("plots/grm_recap_scatter.pdf")

# Genetic relatedness heatmap
np.fill_diagonal(grm_mat, np.nan)
plt.clf()
hm = sns.heatmap(grm_mat, square = True)
hm.get_figure().savefig("plots/grm.pdf")

# Pedigree relatedness heatmap
np.fill_diagonal(prm_mat, np.nan)
plt.clf()
pm = sns.heatmap(prm_mat, square= True)
pm.get_figure().savefig("plots/prm.pdf")

