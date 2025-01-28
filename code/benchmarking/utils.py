import itertools
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator,eigsh

def grm_summary_stat(ts, sample_sets, mode="branch"):
    
    n = len(sample_sets)
    indexes = [(n1, n2) for n1, n2 in itertools.combinations_with_replacement(range(n), 2)]
    K = np.zeros((n,n))
    K[np.triu_indices(n)] = ts.genetic_relatedness(sample_sets, indexes, mode = mode, proportion=False, span_normalise=False)
    K = K + np.triu(K,1).transpose()
    return K

def ts_grm_pca(ts, mode, n_pc=10):
    sample_sets = [[i] for i in ts.samples()]
    grm = ts.genetic_relatedness_matrix(sample_sets, mode=mode)
    eigval, eigvec = eigh(grm, subset_by_index=[grm.shape[0] - n_pc, grm.shape[0] - 1])
    prin_comp = (grm @ eigvec[:, ::-1])
    return prin_comp

def ts_linop_pca(ts, mode, n_pc=10):
    grm_matvec = lambda x: ts.genetic_relatedness_vector(x, mode=mode)
    grm = LinearOperator((ts.num_samples, ts.num_samples), matvec=grm_matvec)
    eigval, eigvec = eigsh(grm, k=10)
    prin_comp = (grm @ eigvec[:, ::-1])
    return prin_comp

