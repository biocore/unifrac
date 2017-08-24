#distutils: language = c++
from libcpp.unordered_map cimport unordered_map
from libcpp cimport bool

cdef extern from "../sucpp/api.hpp" namespace "su":
    struct mat:
        double* condensed_form
        unsigned int n_samples
        unsigned int cf_size
        char** sample_ids

    enum compute_status:
        okay, 
        tree_missing,
        table_missing,
        unknown_method

    compute_status one_off(const char* biom_filename, const char* tree_filename, 
                               const char* unifrac_method, bool variance_adjust, double alpha,
                               unsigned int threads, mat* &result)
    void destroy_mat(mat* &result)
