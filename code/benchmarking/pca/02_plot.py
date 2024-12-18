import numpy as np
import pandas as pd

from plotnine import *

# Aggregate simulations
out_dir = "output/simulations/pca/"
sim_df = pd.read_csv(out_dir + "sim0.csv")
n_sims = 51
for i in range(1, n_sims):
    this_df = pd.read_csv(out_dir + "sim" + str(i) + ".csv")
    sim_df = pd.concat([sim_df, this_df], axis=0)

# Simulated data parameters
default_sample_size = 2**9
default_sequence_length = 2**19
default_mutation_rate = 10**-8

#sim_df = pd.read_csv("output/benchmarking/simulation_results.csv",
#                     names=["rep", "num_samples", "seq_length", "mutation_rate", "method", "time", "mode"])

sim_df["log_mut_rate"] = np.log(sim_df["mutation_rate"])
sim_df["log_time"] = np.log(sim_df["time"])
sim_df["log2_length"] = np.log2(sim_df["seq_length"].astype(int))
sim_df["log2_num_samples"] = np.log2(sim_df["num_samples"].astype(int))

sim_df = sim_df[sim_df["mode"] == "site"]
sim_df = sim_df.drop(["mode"], axis = 1)
### PLOT NUMBER OF SAMPLES ###
plot_dir = "plots/simulations/pca/"

plot_df = sim_df[(sim_df["seq_length"] == default_sequence_length) & (sim_df["mutation_rate"] == default_mutation_rate)]
#plot_df = plot_df.groupby(["method", "mode", "log2_num_samples"]).mean().reset_index()
plot_df = plot_df.groupby(["method", "log2_num_samples"]).mean().reset_index()


plt = (ggplot(plot_df) + 
 aes(x = 'log2_num_samples', y = 'log_time', color = 'method') +
 geom_point() +
 geom_line() 
# facet_grid('mode ~ .')
)

plt.save(plot_dir + "num_samples.png")

### PLOT SEQUENCE LENGTH ###
plot_df = sim_df[(sim_df["num_samples"] == default_sample_size) & (sim_df["mutation_rate"] == default_mutation_rate)]
#plot_df = plot_df.groupby(["method", "mode", "log2_length"]).mean().reset_index()
plot_df = plot_df.groupby(["method", "log2_length"]).mean().reset_index()


plt = (ggplot(plot_df) + 
 aes(x = 'log2_length', y = 'log_time', color = 'method') +
 geom_point() +
 geom_line() 
# facet_grid('mode ~ .')
)

plt.save(plot_dir + "seq_length.png")

### PLOT MUTATION RATE ###

plot_df = sim_df[(sim_df["num_samples"] == default_sample_size) & (sim_df["seq_length"] == default_sequence_length)]
#plot_df = plot_df.groupby(["method", "mode", "log_mut_rate"]).mean().reset_index()
plot_df = plot_df.groupby(["method", "log_mut_rate"]).mean().reset_index()


plt = (ggplot(plot_df) + 
 aes(x = 'log_mut_rate', y = 'log_time', color = 'method') +
 geom_point() +
 geom_line() 
# facet_grid('mode ~ .')
)

plt.save(plot_dir + "mutation_rate.png")

