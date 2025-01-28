import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

prm_file = "output/chr3_prm.csv"
prm_df = pd.read_csv(prm_file)

ped_df = pd.read_csv("data/balsac_pedigree.csv")
ind_dict = dict(zip(ped_df.pedID, ped_df.ind))
new_cols = [str(ind_dict[int(i)]) for i in prm_df.columns.values]
prm_df.columns = new_cols
prm_df = prm_df.set_index(prm_df.columns.values)

# Get selected individuals in region order
selected_df = pd.read_csv("data/balsac_subsample.csv")
order_cols = selected_df.sort_values("proband_region")['ind'].tolist()
order_cols = [str(i) for i in order_cols]

# Get relative information
relative_df = pd.read_csv("data/balsac_relatives.csv")
relative_df = relative_df.loc[relative_df['proband1'].isin(selected_df['ind']) & relative_df['proband2'].isin(selected_df['ind'])]

plot_pairs = (
    relative_df.loc[relative_df['proband_region1'] == relative_df['proband_region2']]
    .groupby(['proband_region1', 'relationship'])
    .head(2)
)
plot_pairs['proband1'] = plot_pairs['proband1'].astype(str)
plot_pairs['proband2'] = plot_pairs['proband2'].astype(str)
plot_pairs['prm'] = plot_pairs.apply(lambda x: prm_df.loc[x['proband2'], x['proband1']], axis=1)

grm_dir = "grms/balsac/chr3/"
full_df = pd.DataFrame()
for type in ["branch", "branch_recap", "site", "site_recap"]:
    plot_df = plot_pairs.copy()
    for sim in range(1,101):
        grm_file = grm_dir + "sim" + str(sim) + "/" + type + ".csv"
        try:
            grm_df = pd.read_csv(grm_file)
            grm_df = grm_df.set_index(grm_df.columns.values)
            plot_df[str(sim)] = plot_pairs.apply(lambda x: grm_df.loc[x['proband2'], x['proband1']], axis=1)
        except:
            pass
    melt_df = pd.melt(plot_df, id_vars=list(plot_pairs.columns.values),
                    var_name = "sim", value_name = "grm")
    melt_df['type'] = type
    full_df = pd.concat([full_df, melt_df])

share = True
if share:
    all_probands = list(set(list(full_df.proband1)) | set(list(full_df.proband2)))
    scramble_map = {j:i for i,j in enumerate(all_probands)}
    full_df['proband1'] = [scramble_map[i] for i in full_df['proband1']]
    full_df['proband2'] = [scramble_map[i] for i in full_df['proband2']]
    full_df.to_csv("output/relatives_relatedness.csv")


plt.clf()
fig, ax1 = plt.subplots(figsize = (8,4))
sns.boxplot(x="prm", y = "grm", hue= "type", data = full_df.loc[full_df['type'].isin(["branch", "site"])])
ax1.set_xticklabels(['{:.2f}'.format(float(t.get_text())) for t in ax1.get_xticklabels()])
plt.xticks(rotation=90)
fig.savefig("plots/boxplot_decap.pdf")

plt.clf()
fig, ax1 = plt.subplots(figsize = (8,4))
sns.boxplot(x="prm", y = "grm", hue= "proband_region1", data = full_df.loc[full_df['type'].isin(["branch_recap"])])
ax1.set_xticklabels(['{:.2f}'.format(float(t.get_text())) for t in ax1.get_xticklabels()])
plt.xticks(rotation=45)
fig.savefig("plots/boxplot_branch_region.pdf")

plt.clf()
fig, ax1 = plt.subplots(figsize = (8,4))
g = sns.FacetGrid(full_df.loc[full_df['type'].isin(["branch_recap"])], col = "generation")
g.map(sns.boxplot, x="prm", y="grm", hue= "proband_region1", data = full_df.loc[full_df['type'].isin(["branch_recap"])])
ax1.set_xticklabels(['{:.2f}'.format(float(t.get_text())) for t in ax1.get_xticklabels()])
plt.xticks(rotation=45)
fig.savefig("plots/boxplot_branch_region_split.pdf")

plt.clf()
fig, ax1 = plt.subplots(figsize = (8,4))
sns.catplot(
    x="prm", y = "grm", 
    kind = "box", col = "generation", 
    hue= "proband_region1", 
    data = full_df.loc[full_df['type'].isin(["branch_recap"])]
    )
#ax1.set_xticklabels(['{:.2f}'.format(float(t.get_text())) for t in ax1.get_xticklabels()])
#plt.xticks(rotation=45)
fig.savefig("plots/boxplot_branch_region.pdf")