import pandas as pd
import numpy as np

# Define the parameter values
methods = ['scikit-allel', 'ts.pca', 'scipy-eigsh', 'scipy-eigh']
modes = ['branch', 'site']
sample_sizes = [2**7, 2**8, 2**9, 2**10, 2**11]
sequence_lengths = [2**15, 2**17, 2**19, 2**21, 2**23]
mutation_rates = [10**-6, 10**-7, 10**-8, 10**-9]

# Create a dictionary to hold the parameter values
parameters = {
    'Method': methods,
    'Mode': modes,
    'Number of samples': sample_sizes,
    'Sequence length': sequence_lengths,
    'Mutation rate (site only)': mutation_rates
}

# Set the default values
default_sample_size = 2**9
default_sequence_length = 2**19
default_mutation_rate = 10**-8

# Create a list of dictionaries for each combination of parameter values
data = []
for method in methods:
    for mode in modes:
        if mode == "branch" and method == "scikit-allel":
             continue
        if mode == 'site':
            if method in ["ts.pca","scipy-eigsh"]:
                 continue
            for mutation_rate in mutation_rates:
                data.append({
                    'Method': method,
                    'Mode': mode,
                    'Num_samples': default_sample_size,
                    'Sequence_length': default_sequence_length,
                    'Mutation_rate': mutation_rate
                })
        for sample_size in sample_sizes:
                data.append({
                    'Method': method,
                    'Mode': mode,
                    'Num_samples': sample_size,
                    'Sequence_length': default_sequence_length,
                    'Mutation_rate': default_mutation_rate
                })
        for sequence_length in sequence_lengths:
            data.append({
                'Method': method,
                'Mode': mode,
                'Num_samples': default_sample_size,
                'Sequence_length': sequence_length,
                'Mutation_rate': default_mutation_rate
            })

# Create the pandas data frame
df = pd.DataFrame(data)
df = df.drop_duplicates()

# Save the data frame
df.to_csv('data/pca_simulation_parameters.csv', index=False)
