import skbio
import numpy as np
cimport numpy as np

def ssu(const char* biom_filename, const char* tree_filename, 
        const char* unifrac_method, bool variance_adjust, double alpha,
        unsigned int threads):
    cdef:
        mat *result;
        compute_status status;
        np.double_t[:, :] view
        np.ndarray[np.double_t, ndim=2] numpy_arr
        double **cf
        int i
        list ids

    status = one_off(biom_filename, tree_filename, unifrac_method, variance_adjust, alpha, threads, result)

    ids = []
    numpy_arr = np.zeros((result.n_samples, result.n_samples), dtype=np.double)
    for i in range(result.n_samples):
        numpy_arr[i, :] = <np.double_t[:result.n_samples]> result.condensed_form[i]
        ids.append(result.sample_ids[i].decode('utf-8'))

    return skbio.DistanceMatrix(numpy_arr, ids)
    
