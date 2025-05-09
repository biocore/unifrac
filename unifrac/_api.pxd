#cython: language_level=3

from libc.stdint cimport uint32_t, uint8_t
from libcpp cimport bool


cdef extern from "status_enum.hpp":
    enum compute_status:
        okay, 
        tree_missing,
        table_missing,
        table_empty,
        unknown_method,
        table_and_tree_do_not_overlap,
        output_error,
        invalid_method,
        grouping_missing


cdef extern from "api.hpp":
    struct mat:
        double* condensed_form
        unsigned int n_samples
        unsigned int cf_size
        char** sample_ids

    struct mat_full_fp64:
        uint32_t n_samples
        uint32_t flags 
        double* matrix
        char** sample_ids
    
    struct mat_full_fp32:
        uint32_t n_samples
        uint32_t flags 
        float* matrix
        char** sample_ids

    struct results_vec:
        unsigned int n_samples
        double* values
        char** sample_ids

    struct support_biom:
        char** obs_ids
        char** sample_ids
        uint32_t* indices
        uint32_t* indptr
        double* data
        int n_obs
        int n_samples
        int nnz

    struct support_bptree:
        bool* structure
        double* lengths
        char** names
        int n_parens

    compute_status one_off(const char* biom_filename, const char* tree_filename,
                           const char* unifrac_method, bool variance_adjust, double alpha,
                           bool bypass_tips, unsigned int n_substeps, mat** result)
    
    compute_status one_off_matrix(const char* biom_filename, const char* tree_filename,
                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                  bool bypass_tips, unsigned int n_substeps, const char *mmap_dir,
                                  mat_full_fp64** result)

    compute_status one_off_matrix_fp32(const char* biom_filename, const char* tree_filename,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, unsigned int n_substeps, const char *mmap_dir,
                                       mat_full_fp32** result)

    compute_status one_off_inmem(const support_biom *table, const support_bptree *tree,
                                 const char* unifrac_method, bool variance_adjust, double alpha,
                                 bool bypass_tips, unsigned int n_substeps, mat_full_fp64** result)

    compute_status one_off_inmem_fp32(const support_biom *table, const support_bptree *tree,
                                      const char* unifrac_method, bool variance_adjust, double alpha,
                                      bool bypass_tips, unsigned int n_substeps, mat_full_fp32** result)

    compute_status one_dense_pair_v2(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
                                       const support_bptree *tree,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, double* result);


    compute_status faith_pd_one_off(const char* biom_filename, const char* tree_filename,
                                    results_vec** result)

    void destroy_mat(mat** result)
    void destroy_mat_full_fp32(mat_full_fp32** result)
    void destroy_mat_full_fp64(mat_full_fp64** result)
    void destroy_results_vec(results_vec** result)

    void ssu_set_random_seed(unsigned int new_seed)

    compute_status unifrac_to_file_v2(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                     const char* unifrac_method, bool variance_adjust, double alpha,
                                     bool bypass_tips, unsigned int n_substeps, const char* format,
                                     unsigned int subsample_depth, bool subsample_with_replacement,
                                     unsigned int pcoa_dims,
                                     unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                     const char *mmap_dir)

    # obsolete, will be removed in the future
    compute_status unifrac_to_file(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                     const char* unifrac_method, bool variance_adjust, double alpha,
                                     bool bypass_tips, unsigned int n_substeps, const char* format,
                                     unsigned int pcoa_dims, const char *mmap_dir)

    compute_status unifrac_multi_to_file_v2(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                            const char* unifrac_method, bool variance_adjust, double alpha,
                                            bool bypass_tips, unsigned int n_substeps, const char* format,
                                            unsigned int n_subsamples, unsigned int subsample_depth, bool subsample_with_replacement,
                                            unsigned int pcoa_dims,
                                            unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                            const char *mmap_dir)


