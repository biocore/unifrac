#cython: language_level=3
#distutils: language = c++

from libcpp cimport bool, string
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport uint32_t


cdef extern from "biom.hpp" namespace "su":
    cdef cppclass biom:
        biom(const vector[string] &obs_ids, 
             const vector[string] &samp_ids, 
             const vector[uint32_t] &index,
             const vector[uint32_t] &indptr,
             const vector[double] &data)
        biom()

cdef extern from "tree.hpp" namespace "su":
    cdef cppclass BPTree:
        BPTree(vector[bool] input_structure,
               vector[double] input_lengths,
               vector[string] input_names)
        BPTree()

cdef extern from "api.hpp":
    struct mat:
        double* condensed_form
        unsigned int n_samples
        unsigned int cf_size
        char** sample_ids

    struct results_vec:
        unsigned int n_samples
        double* values
        char** sample_ids

    enum compute_status:
        okay, 
        tree_missing,
        table_missing,
        table_empty,
        unknown_method,
        table_and_tree_do_not_overlap,
        output_error

    compute_status one_off(const char* biom_filename, const char* tree_filename, 
                               const char* unifrac_method, bool variance_adjust, double alpha,
                               bool bypass_tips, unsigned int threads, mat** result)
    
    compute_status one_off_inmem(biom *table, BPTree *tree, 
                                 const char* unifrac_method, bool variance_adjust, double alpha,
                                 bool bypass_tips, unsigned int threads, mat** result)

    compute_status faith_pd_one_off(const char* biom_filename, const char* tree_filename,
                                    results_vec** result)

    void destroy_mat(mat** result)

    void destroy_results_vec(results_vec** result)

    compute_status unifrac_to_file(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                     const char* unifrac_method, bool variance_adjust, double alpha,
                                     bool bypass_tips, unsigned int threads, const char* format,
                                     unsigned int pcoa_dims, const char *mmap_dir)

