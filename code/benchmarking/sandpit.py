import sys

sys.path.insert(0, "/home/brieuc/Documents/Projects/tskit/python")
print(sys.path)
import tskit
print(tskit.__version__)
print(tskit.__file__)

import msprime

num_samples = 2**5
seq_length = 2**11
mutation_rate = 10**-8

seed = 42                   # Random seed
recombination_rate = 1e-8   # Recombination rate per base
Ne = 1e4                    # Effective population size

num_simulations = 10        # Number of simulations to perform

import itertools
import numpy as np
from egrm import varGRM_C

ts = msprime.simulate(num_samples,
                      length=seq_length,
                      recombination_rate=recombination_rate,
                      Ne=Ne,
                      mutation_rate=mutation_rate,
                      random_seed=seed)

sample_sets = [[i] for i in ts.samples()]

def grm_summary_stat(ts, sample_sets, mode="branch"):
    
    n = len(sample_sets)
    indexes = [(n1, n2) for n1, n2 in itertools.combinations_with_replacement(range(n), 2)]
    K = np.zeros((n,n))
    K[np.triu_indices(n)] = ts.genetic_relatedness(sample_sets, indexes, mode = mode, proportion=False, span_normalise=False)
    K = K + np.triu(K,1).transpose()
    return K


egrm = varGRM_C(ts, var=False)
divmat = ts.genetic_relatedness_matrix(sample_sets)
ts_grm = grm_summary_stat(ts)

num_samples = 2**5
seq_length = 2**11
mutation_rate = 10**-8

seed = 42                   # Random seed
recombination_rate = 1e-8   # Recombination rate per base
Ne = 1e4                    # Effective population size

num_simulations = 10   

ts_sim = msprime.simulate(num_samples,
                      length=seq_length,
                      recombination_rate=recombination_rate,
                      Ne=Ne,
                      mutation_rate=mutation_rate,
                      random_seed=seed)
sample_sets = [[i] for i in ts_sim.samples()]
divmat = ts_sim.genetic_relatedness_matrix(sample_sets, mode = "branch", span_normalise=False)