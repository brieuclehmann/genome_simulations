
import tskit
import msprime
import sys
import os

sys.path.insert(0, "")

### Import packages ###
import pandas as pd
import timeit

from egrm import varGRM_C
from utils import grm_summary_stat

### Set / get simulation parameters ### 
rep = int(sys.argv[1]) - 1

param_df = pd.read_csv('data/grm_simulation_parameters.csv')
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
    sample_sets = [[i] for i in ts.samples()]
    # Measure the time for GRM computation
    if method == "eGRM":
        grm_time = timeit.timeit(lambda: varGRM_C(ts, var=False), number=1)
    elif method == "ts.summary_stat":
        grm_time = timeit.timeit(lambda: grm_summary_stat(ts, sample_sets, mode=mode), number=1)
    elif method == "ts.genetic_relatedness_matrix":
        grm_time = timeit.timeit(lambda: ts.genetic_relatedness_matrix(sample_sets, mode=mode), number=1)
    # Concatenate results
    sim_df = pd.concat([sim_df, pd.DataFrame([{"rep": i,
                                               "num_samples": num_samples,
                                               "seq_length": seq_length,
                                               "mutation_rate": mutation_rate,
                                               "method": method,
                                               "mode": mode,
                                               "time": grm_time}])])

# Save results
outdir = "output/benchmarking/grm"
if not os.path.exists(outdir):
    os.makedirs(outdir)

outfile = outdir + "/sim" + str(rep) + ".csv"
sim_df.to_csv(outfile, index=False)
