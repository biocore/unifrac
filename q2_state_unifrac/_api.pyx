import skbio
import numpy as np
cimport numpy as np

def ssu(str biom_filename, str tree_filename, 
        str unifrac_method, bool variance_adjust, double alpha,
        unsigned int threads):
    cdef:
        mat *result;
        compute_status status;
        np.ndarray[np.double_t, ndim=1] numpy_arr
        double *cf
        int i
        bytes biom_py_bytes
        bytes tree_py_bytes
        bytes met_py_bytes
        char* biom_c_string 
        char* tree_c_string 
        char* met_c_string 
        list ids

    biom_py_bytes = biom_filename.encode()
    tree_py_bytes = tree_filename.encode()
    met_py_bytes = unifrac_method.encode()
    biom_c_string = biom_py_bytes
    tree_c_string = tree_py_bytes
    met_c_string = met_py_bytes

    status = one_off(biom_c_string, 
                     tree_c_string, 
                     met_c_string, 
                     variance_adjust, 
                     alpha, 
                     threads, 
                     &result)

    if status != okay:
        if status == tree_missing:
            raise ValueError("Tree file not found.")
        if status == table_missing:
            raise ValueError("Table file not found.")
        if status == unknown_method:
            raise ValueError("Unknown method.")

    ids = []
    numpy_arr = np.zeros(result.cf_size, dtype=np.double)
    numpy_arr[:] = <np.double_t[:result.cf_size]> result.condensed_form
    
    for i in range(result.n_samples):
        ids.append(result.sample_ids[i].decode('utf-8'))

    destroy_mat(&result)

    return skbio.DistanceMatrix(numpy_arr, ids)
    
