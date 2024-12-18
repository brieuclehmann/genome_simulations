#TODO: Check site implementation for eGRM
#TODO: Check call to ts.genetic_relatedness

import sys
import os
import allel

sys.path.insert(0, "scripts/benchmarking")
import utils

sys.path.insert(0, "../tskit/python")
import tskit
import msprime
print(tskit.__version__)

### Import packages ###
import pandas as pd
import timeit

### Set / get simulation parameters ### 
rep = int(sys.argv[1]) - 1

param_df = pd.read_csv('data/pca_simulation_parameters.csv')
this_sim = param_df.loc[rep]

num_samples = int(this_sim["Num_samples"])          # Number of samples (individuals)
seq_length =  int(this_sim["Sequence_length"])      # Length of the genomic sequence
mutation_rate = float(this_sim["Mutation_rate"])    # Mutation rate per site
method = this_sim["Method"]                         # Method to use for computing the GRM
mode = this_sim["Mode"]                             # Mode (site v branch) to use for computing the GRM

seed = 42                   # Random seed
recombination_rate = 1e-8   # Recombination rate per base
Ne = 1e4                    # Effective population size

num_simulations = 10        # Number of simulations to perform

### Perform the simulation study ###
sim_df = pd.DataFrame(columns=["rep", "num_samples", "seq_length", "mutation_rate", "method", "time"])

for i in range(num_simulations):
    # Simulate a tree sequence with mutations
    ts = msprime.simulate(num_samples,
                            length=seq_length,
                            recombination_rate=recombination_rate,
                            Ne=Ne,
                            mutation_rate=mutation_rate,
                            random_seed=seed+i)
    # Measure the time for GRM computation
    if method == "scikit-allel":
        gn = ts.genotype_matrix()
        pca_time = timeit.timeit(lambda: allel.randomized_pca(gn, iterated_power=5, ploidy=1), number=1)
    elif method == "ts.pca":
        pca_time = timeit.timeit(lambda: ts.pca(num_components=10, mode=mode), number=1)
    elif method == "scipy-eigh":
        pca_time = timeit.timeit(lambda: utils.ts_grm_pca(ts, mode=mode), number=1)
    elif method == "scipy-eigsh":
        pca_time = timeit.timeit(lambda: utils.ts_linop_pca(ts, mode=mode), number=1)
    # Concatenate results
    sim_df = pd.concat([sim_df, pd.DataFrame([{"rep": i,
                                               "num_samples": num_samples,
                                               "seq_length": seq_length,
                                               "mutation_rate": mutation_rate,
                                               "method": method,
                                               "mode": mode,
                                               "time": pca_time}])])

# Save results
outdir = "output/simulations/pca"
if not os.path.exists(outdir):
    os.makedirs(outdir)

outfile = outdir + "/sim" + str(rep) + ".csv"
sim_df.to_csv(outfile, index=False)
