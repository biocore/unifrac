#cimport libc.stdlib
#from cpython cimport Py_INCREF, Py_DECREF, bool
cimport numpy as np
np.import_array()

def ssu(const char* biom_filename, const char* tree_filename, 
        const char* unifrac_method, bool variance_adjust, double alpha,
        unsigned int threads):
    cdef:
        mat* result;
        compute_status status;

    status = one_off(biom_filename, tree_filename, unifrac_method, variance_adjust, alpha, threads, result)
    print(status)
    
