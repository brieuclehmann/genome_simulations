import time

def compute_grm(ts, grm_ind, mode = "branch"):

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

    start = time.time()
    grm = ts_sub.genetic_relatedness_matrix(sample_sets, mode = mode)
    end = time.time()
    print(end - start)

    return grm

def simplify_ts(ts, grm_ind):
    grm_sample_sets = []
    for i in grm_ind:
        grm_sample_sets.append([i.nodes[0]])
        grm_sample_sets.append([i.nodes[1]])
    ts_sub = ts.simplify([i[0] for i in grm_sample_sets], filter_individuals=False)
    return ts_sub