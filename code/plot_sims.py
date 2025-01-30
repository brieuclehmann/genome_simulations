import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

prm_file = "output/chr3_prm.csv"
prm_df = pd.read_csv(prm_file)

prm_df['proband1'] = prm_df.columns.values
prm_long = pd.melt(prm_df, id_vars='proband1', var_name='proband2', value_name='prm')
prm_long = prm_long[prm_long['proband1'] <= prm_long['proband2']]

# Get selected individuals in region order
selected_df = pd.read_csv("data/balsac_subsample.csv")
order_cols = selected_df.sort_values("proband_region")['ind'].tolist()
order_cols = [str(i) for i in order_cols]

# Get relative information
relative_df = pd.read_csv("data/balsac_relatives.csv")
relative_df1 = relative_df.loc[relative_df['proband1'].isin(selected_df['ind']) & relative_df['proband2'].isin(selected_df['ind'])]
relative_df1 = relative_df1[['proband1', 'proband2', 'generation', 'relationship']]
relative_df2 = relative_df1.rename(columns={'proband1':'proband2','proband2':'proband1'})
relative_df = pd.concat([relative_df1,relative_df2])
relative_df = relative_df[relative_df['proband1'] <= relative_df['proband2']]

grm_dir = "grms/balsac/chr3/"
full_df = pd.DataFrame()
for type in ["branch_recap"]:
    for sim in range(1,101):
        grm_file = grm_dir + "sim" + str(sim) + "/" + type + ".csv"
        try:
            grm_df = pd.read_csv(grm_file)
            grm_df['proband1'] = grm_df.columns.values
            # Reset index to make row names a column
            melted = pd.melt(grm_df, id_vars='proband1', var_name='proband2', value_name='grm')
            melted = melted[melted['proband1'] <= melted['proband2']]
            melted['sim'] = sim
            melted['type'] = type
            full_df = pd.concat([full_df, melted])
        except:
            pass

full_df = full_df.merge(prm_long, on = ["proband1", "proband2"])
full_df = full_df.astype({"proband1":"int64", "proband2":"int64"})
full_df = full_df.merge(relative_df, on = ["proband1", "proband2"], how = "left")

full_df.to_csv("output/relatives_relatedness.csv")

share = True
if share:
    meta_df = pd.read_csv("data/balsac_subsample_meta.csv")
    #all_probands = list(set(list(full_df.proband1)) | set(list(full_df.proband2)))
    #scramble_map = {j:i for i,j in enumerate(all_probands)}
    share_df = full_df.merge(meta_df, left_on = "proband1", right_on = "proband")
    share_df = share_df.drop(columns = ["proband", "proband1"])
    share_df = share_df.rename(
        columns = {"proband_region":"proband_region1",
                    "noid":"proband1",
                    "average_depth":"average_depth1"}
                    )
    share_df = share_df.merge(meta_df, left_on = "proband2", right_on = "proband")
    share_df = share_df.drop(columns = ["proband", "proband2"])
    share_df = share_df.rename(
        columns = {"proband_region":"proband_region2",
                    "noid":"proband2",
                    "average_depth":"average_depth2"}
                    )
    # full_df['proband1'] = [scramble_map[i] for i in full_df['proband1']]
    # full_df['proband2'] = [scramble_map[i] for i in full_df['proband2']]
    share_df.to_csv("output/relatives_relatedness_noid.csv", index = False)


# plt.clf()
# fig, ax1 = plt.subplots(figsize = (8,4))
# sns.boxplot(x="prm", y = "grm", hue= "type", data = full_df.loc[full_df['type'].isin(["branch", "site"])])
# ax1.set_xticklabels(['{:.2f}'.format(float(t.get_text())) for t in ax1.get_xticklabels()])
# plt.xticks(rotation=90)
# fig.savefig("plots/boxplot_decap.pdf")

# plt.clf()
# fig, ax1 = plt.subplots(figsize = (8,4))
# sns.boxplot(x="prm", y = "grm", hue= "proband_region1", data = full_df.loc[full_df['type'].isin(["branch_recap"])])
# ax1.set_xticklabels(['{:.2f}'.format(float(t.get_text())) for t in ax1.get_xticklabels()])
# plt.xticks(rotation=45)
# fig.savefig("plots/boxplot_branch_region.pdf")

# plt.clf()
# fig, ax1 = plt.subplots(figsize = (8,4))
# g = sns.FacetGrid(full_df.loc[full_df['type'].isin(["branch_recap"])], col = "generation")
# g.map(sns.boxplot, x="prm", y="grm", hue= "proband_region1", data = full_df.loc[full_df['type'].isin(["branch_recap"])])
# ax1.set_xticklabels(['{:.2f}'.format(float(t.get_text())) for t in ax1.get_xticklabels()])
# plt.xticks(rotation=45)
# fig.savefig("plots/boxplot_branch_region_split.pdf")

# plt.clf()
# fig, ax1 = plt.subplots(figsize = (8,4))
# sns.catplot(
#     x="prm", y = "grm", 
#     kind = "box", col = "generation", 
#     hue= "proband_region1", 
#     data = full_df.loc[full_df['type'].isin(["branch_recap"])]
#     )
# #ax1.set_xticklabels(['{:.2f}'.format(float(t.get_text())) for t in ax1.get_xticklabels()])
# #plt.xticks(rotation=45)
# fig.savefig("plots/boxplot_branch_region.pdf")