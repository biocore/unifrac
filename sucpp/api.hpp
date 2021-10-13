#include "task_parameters.hpp"

#ifdef __cplusplus
#include <vector>
#define EXTERN extern "C"


#else
#include <stdbool.h>
#define EXTERN
#endif

#define PARTIAL_MAGIC "SSU-PARTIAL-01"
#define PARTIAL_MAGIC_V2 0x088ABA02


typedef enum compute_status {okay=0, tree_missing, table_missing, table_empty, unknown_method, table_and_tree_do_not_overlap, output_error} ComputeStatus;
typedef enum io_status {read_okay=0, write_okay, open_error, read_error, magic_incompatible, bad_header, unexpected_end, write_error} IOStatus;
typedef enum merge_status {merge_okay=0, incomplete_stripe_set, sample_id_consistency, square_mismatch, partials_mismatch, stripes_overlap} MergeStatus;

/* a result matrix
 *
 * n_samples <uint> the number of samples.
 * cf_size <uint> the size of the condensed form.
 * is_upper_triangle <bool> if true, indicates condensed_form represents a square
 *      matrix, and only the upper triangle is contained. if false,
 *      condensed_form represents the lower triangle of a matrix.
 * condensed_form <double*> the matrix values of length cf_size.
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct mat {
    unsigned int n_samples;
    unsigned int cf_size;
    bool is_upper_triangle;
    double* condensed_form;
    char** sample_ids;
} mat_t;

/* a result matrix, full, fp64
 *
 * n_samples <uint> the number of samples.
 * matrix <double*> the matrix values, n_sample**2 size
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct mat_full_fp64 {
    uint32_t n_samples;
    uint32_t flags; //opaque, 0 for default behavior
    double* matrix;
    char** sample_ids;
} mat_full_fp64_t;

/* a result matrix, full, fp32
 *
 * n_samples <uint> the number of samples.
 * matrix <float*> the matrix values, n_sample**2 size
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct mat_full_fp32 {
    uint32_t n_samples;
    uint32_t flags; //opaque, 0 for default behavior
    float* matrix;
    char** sample_ids;
} mat_full_fp32_t;



/* a result vector
 *
 * n_samples <uint> the number of samples.
 * values <double*> the score values of length n_samples.
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct results_vec{
    unsigned int n_samples;
    double* values;
    char** sample_ids;
} r_vec;

/* a partial result containing stripe data
 *
 * n_samples <uint> the number of samples.
 * sample_ids <char**> the sample IDs of length n_samples.
 * stripes <double**> the stripe data of dimension (stripe_stop - stripe_start, n_samples)
 * stripe_start <uint> the logical starting stripe in the final matrix.
 * stripe_stop <uint> the logical stopping stripe in the final matrix.
 * stripe_total <uint> the total number of stripes present in the final matrix.
 * is_upper_triangle <bool> whether the stripes correspond to the upper triangle of the resulting matrix.
 *      This is useful for asymmetric unifrac metrics.
 */
typedef struct partial_mat {
    uint32_t n_samples;
    char** sample_ids;
    double** stripes;
    uint32_t stripe_start;
    uint32_t stripe_stop;
    uint32_t stripe_total;
    bool is_upper_triangle;
} partial_mat_t;

/* a partial resuly, can be populated dynamically
 *
 * n_samples <uint> the number of samples.
 * sample_ids <char**> the sample IDs of length n_samples.
 * offsets <uint64_t*> offsets to the stripes in the file; 0 means unknown
 * stripes <double**> the stripe data of dimension (stripe_stop - stripe_start, n_samples)
 * stripe_start <uint> the logical starting stripe in the final matrix.
 * stripe_stop <uint> the logical stopping stripe in the final matrix.
 * stripe_total <uint> the total number of stripes present in the final matrix.
 * is_upper_triangle <bool> whether the stripes correspond to the upper triangle of the resulting matrix.
 *      This is useful for asymmetric unifrac metrics.
 * filename <char*> Name of the file from which to read
 */
typedef struct partial_dyn_mat {
    uint32_t n_samples;
    char** sample_ids;
    uint64_t* offsets;
    double** stripes;
    uint32_t stripe_start;
    uint32_t stripe_stop;
    uint32_t stripe_total;
    bool is_upper_triangle;
    char* filename;
} partial_dyn_mat_t;



void destroy_mat(mat_t** result);
void destroy_mat_full_fp64(mat_full_fp64_t** result);
void destroy_mat_full_fp32(mat_full_fp32_t** result);
void destroy_partial_mat(partial_mat_t** result);
void destroy_partial_dyn_mat(partial_dyn_mat_t** result);
void destroy_results_vec(r_vec** result);

/* Compute UniFrac - condensed form
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * threads <uint> the number of threads to use.
 * result <mat_t**> the resulting distance matrix in condensed form, this is initialized within the method so using **
 *
 * one_off returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_off(const char* biom_filename, const char* tree_filename,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, unsigned int threads, mat_t** result);

/* Compute UniFrac - matrix form
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * threads <uint> the number of threads/blocks to use.
 * mmap_dir <const char*> If not NULL, area to use for temp memory storage
 * result <mat_full_fp64_t**> the resulting distance matrix in matrix form, this is initialized within the method so using **
 *
 * one_off_matrix returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_off_matrix(const char* biom_filename, const char* tree_filename,
                                    const char* unifrac_method, bool variance_adjust, double alpha,
                                    bool bypass_tips, unsigned int nthreads,
                                    const char *mmap_dir,
                                    mat_full_fp64_t** result);

/* Compute UniFrac - matrix form, fp32 variant
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * threads <uint> the number of threads/blocks to use.
 * mmap_dir <const char*> If not NULL, area to use for temp memory storage
 * result <mat_full_fp32_t**> the resulting distance matrix in matrix form, this is initialized within the method so using **
 *
 * one_off_matrix_fp32 returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_off_matrix_fp32(const char* biom_filename, const char* tree_filename,
                                         const char* unifrac_method, bool variance_adjust, double alpha,
                                         bool bypass_tips, unsigned int nthreads,
                                         const char *mmap_dir,
                                         mat_full_fp32_t** result);


/* compute Faith PD
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * result <r_vec**> the resulting vector of computed Faith PD values
 *
 * faith_pd_one_off returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus faith_pd_one_off(const char* biom_filename, const char* tree_filename,
                                      r_vec** result);

/* Compute UniFrac and save to file
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * out_filename <const char*> the filename of the output file.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * threads <uint> the number of threads to use.
 * format <const char*> output format to use.
 * pcoa_dims <uint> if not 0, number of dimensions to use or PCoA
 * mmap_dir <const char*> if not empty, temp dir to use for disk-based memory 
 *
 * unifrac_to_file returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 * output_error   : failed to properly write the output file
 */
EXTERN ComputeStatus unifrac_to_file(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                     const char* unifrac_method, bool variance_adjust, double alpha,
                                     bool bypass_tips, unsigned int threads, const char* format,
                                     unsigned int pcoa_dims, const char *mmap_dir);

/* Write a matrix object
 *
 * filename <const char*> the file to write into
 * result <mat_t*> the results object
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_mat(const char* filename, mat_t* result);

/* Write a matrix object using hdf5 format
 *
 * filename <const char*> the file to write into
 * result <mat_t*> the results object
 * pcoa_dims <uint> PCoAdimensions to compute, if >0
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_mat_hdf5(const char* filename, mat_t* result, unsigned int pcoa_dims);

/* Write a matrix object using hdf5 format, using fp32 precision
 *
 * filename <const char*> the file to write into
 * result <mat_t*> the results object
 * pcoa_dims <uint> PCoAdimensions to compute, if >0
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_mat_hdf5_fp32(const char* filename, mat_t* result, unsigned int pcoa_dims);

/* Write a matrix object
 *
 * filename <const char*> the file to write into
 * result <mat_full_t*> the results object
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_mat_from_matrix(const char* filename, mat_full_fp64_t* result);


/* Write a matrix object from buffer using hdf5 format
 *
 * filename <const char*> the file to write into
 * result <mat_full_t*> the results object
 * pcoa_dims <uint> PCoAdimensions to compute, if >0
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_mat_from_matrix_hdf5(const char* filename, mat_full_fp64_t* result, unsigned int pcoa_dims);

/* Write a matrix object from buffer using hdf5 format, using fp32 precision
 *
 * filename <const char*> the file to write into
 * result <mat_full_fp32_t*> the results object
 * pcoa_dims <uint> PCoAdimensions to compute, if >0
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_mat_from_matrix_hdf5_fp32(const char* filename, mat_full_fp32_t* result, unsigned int pcoa_dims);

/* Write a series
 *
 * filename <const char*> the file to write into
 * result <r_vec*> the results object
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_vec(const char* filename, r_vec* result);

/* Read a matrix object
 *
 * filename <const char*> the file to write into
 * result <mat_t**> the results object
 *
 * The following error codes are returned:
 *
 * read_okay : no problems
 * open_error : could not open the file
 * magic_incompatible : format magic not found or incompatible
 * unexpected_end     : format end not found in expected location
 */
//EXTERN IOStatus read_mat(const char* filename, mat_t** result);

/* Compute a subset of a UniFrac distance matrix
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * threads <uint> the number of threads to use.
 * stripe_start <uint> the starting stripe to compute
 * stripe_stop <uint> the last stripe to compute
 * dm_stripes <vector of double*> the unique branch length stripes. This is expected to be
 *      uninitialized, and is an output parameter.
 * dm_stripes_total <vector of double*> the total branch length stripes. This is expected to be
 *      uninitialized, and is an output parameter.
 * result <partial_mat_t**> the resulting distance matrix in condensed form, this is initialized within the method so using **
 *
 * partial returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 */

EXTERN ComputeStatus partial(const char* biom_filename, const char* tree_filename,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, unsigned int threads, unsigned int stripe_start,
                             unsigned int stripe_stop, partial_mat_t** result);

/* Write a partial matrix object
 *
 * filename <const char*> the file to write into
 * result <partial_mat_t*> the partial results object
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 * open_error : could not open the file
 *
 * The structure of the binary output file is as follows. Newlines added for clarity, but are not stored.
 * The file has logical blocks, but are not explicitly denoted in the format. These logical blocks are
 * just used to improve readability here, and are denoted by ### marks.
 *
 * ### HEADER ###
 * <MAGIC_LEN>          : uint16_t, the length of the magic
 * <MAGIC>              : char, e.g., SSU-PARTIAL-01
 * <N_SAMPLES>          : uint32_t, the number of samples
 * <N_STRIPES>          : uint32_t, the number of stripes represented in this file
 * <STRIPE_START>       : uint32_t, the starting stripe number
 * <STRIPES_TOTAL>      : uint32_t, the total number of stripes in the full matrix
 * <IS_UPPER_TRIANGLE>  : uint8_t, zero is false, nonzero is true
 *
 * ### SAMPLE IDS ###
 * <LEN>                : uint16_t, the length of the next sample ID
 * <SAMPLE_ID[0]>       : LEN bytes, char
 * ...                  : ... repeated <LEN><SAMPLE_ID[i]>
 * <LEN>                : uint16_t, the length of the next sample ID
 * <SAMPLE_ID[N]>       : LEN bytes, char
 *
 * ### STRIPE VALUES; SS -> STRIPE_START, NS -> N_STRIPES
 * <STRIPE[SS][0]>      : double, the first value in the 0th stripe
 * ...                  : ... repeated for N_SAMPLES values
 * <STRIPE[SS][N]>      : double, the last value in the 0th stripe
 * <STRIPE[SS + NS][0]> : double, the first value in the Kth stripe
 * ...                  : ... repeated for N_SAMPLES values
 * <STRIPE[SS + NS][N]> : double, the last value in the Kth stripe
 *
 * ### FOOTER ###
 * <MAGIC>              : char, e.g., SSU-PARTIAL-01, same as starting magic
 */
EXTERN IOStatus write_partial(const char* filename, const partial_mat_t* result);

/* Read a partial matrix object
 *
 * filename <const char*> the file to write into
 * result <partial_mat_t**> the partial results object, output parameter
 *
 * The following error codes are returned:
 *
 * read_okay          : no problems
 * open_error         : could not open the file
 * magic_incompatible : format magic not found or incompatible
 * bad_header         : header seems malformed
 * unexpected_end     : format end not found in expected location
 */
EXTERN IOStatus read_partial(const char* filename, partial_mat_t** result);

/* Read a partial matrix object header
 *
 * filename <const char*> the file to write into
 * result <partial_dyn_mat_t**> the partial results object, output parameter
 *
 * The following error codes are returned:
 *
 * read_okay          : no problems
 * open_error         : could not open the file
 * magic_incompatible : format magic not found or incompatible
 * bad_header         : header seems malformed
 * unexpected_end     : format end not found in expected location
 */
EXTERN IOStatus read_partial_header(const char* input_filename, partial_dyn_mat_t** result_out);

/* Read a stripe of a partial matrix
 *
 * filename <const char*> the file to write into
 * result <partial_dyn_mat_t*> the partial results object
 * stripe_idx <uint> relative stripe number 
 *
 * The following error codes are returned:
 *
 * read_okay          : no problems
 * open_error         : could not open the file
 * magic_incompatible : format magic not found or incompatible
 * bad_header         : header seems malformed
 * unexpected_end     : format end not found in expected location
 */
EXTERN IOStatus read_partial_one_stripe(partial_dyn_mat_t* result, uint32_t stripe_idx);

/*
 * Description TBD
 */
EXTERN MergeStatus validate_partial(const partial_dyn_mat_t* const * partial_mats, int n_partials);

/* Merge partial results
 *
 * results <partial_mat_t**> an array of partial_mat_t*, the buffers will be destroyed in the process
 * n_partials <int> number of partial mats
 * merged <mat_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
EXTERN MergeStatus merge_partial(partial_mat_t** partial_mats, int n_partials, unsigned int nthreads, mat_t** result);

/* Merge partial results
 *
 * partial_mats <partial_dyn_mat_t**> an array of partial_dyn_mat_t*
 * n_partials <int> number of partial mats
 * result <mat_full_fp64_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
MergeStatus merge_partial_to_matrix(partial_dyn_mat_t* * partial_mats, int n_partials, mat_full_fp64_t** result);

/* Merge partial results
 *
 * partial_mats <partial_dyn_mat_t**> an array of partial_dyn_mat_t*
 * n_partials <int> number of partial mats
 * result <mat_full_fp32_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *             
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
MergeStatus merge_partial_to_matrix_fp32(partial_dyn_mat_t* * partial_mats, int n_partials, mat_full_fp32_t** result);


/* Merge partial results
 *
 * partial_mats <partial_dyn_mat_t**> an array of partial_dyn_mat_t*
 * n_partials <int> number of partial mats
 * mmap_dir <const char *> Where to host the mmap file
 * result <mat_full_fp64_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
MergeStatus merge_partial_to_mmap_matrix(partial_dyn_mat_t* * partial_mats, int n_partials, const char *mmap_dir, mat_full_fp64_t** result);

/* Merge partial results
 *
 * partial_mats <partial_dyn_mat_t**> an array of partial_dyn_mat_t*
 * n_partials <int> number of partial mats
 * mmap_dir <const char *> Where to host the mmap file
 * result <mat_full_fp32_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
MergeStatus merge_partial_to_mmap_matrix_fp32(partial_dyn_mat_t* * partial_mats, int n_partials, const char *mmap_dir, mat_full_fp32_t** result);


// Find eigen values and vectors
// Based on N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
//     Original Paper: https://arxiv.org/abs/1007.5510
// centered == n x n, must be symmetric, Note: will be used in-place as temp buffer
void find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, double * centered, double **eigenvalues, double **eigenvectors);
void find_eigens_fast_p32(const uint32_t n_samples, const uint32_t n_dims,float * centered, float **eigenvalues, float **eigenvectors);

// Perform Principal Coordinate Analysis
// mat       - in, result of unifrac compute
// n_samples - in, size of the matrix (n x n)
// n_dims    - in, Dimensions to reduce the distance matrix to. This number determines how many eigenvectors and eigenvalues will be returned.
// eigenvalues - out, alocated buffer of size n_dims
// samples     - out, alocated buffer of size n_dims x n_samples
// proportion_explained - out, allocated buffer of size n_dims
void pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, double **eigenvalues, double **samples, double **proportion_explained);
void pcoa_fp32(const float * mat, const uint32_t n_samples, const uint32_t n_dims, float * *eigenvalues, float * *samples, float * *proportion_explained);
void pcoa_mixed(const double * mat, const uint32_t n_samples, const uint32_t n_dims, float * *eigenvalues, float * *samples, float * *proportion_explained);


#ifdef __cplusplus
// TODO: only needed for testing, should be encased in a macro
void set_tasks(std::vector<su::task_parameters> &tasks,
               double alpha,
               unsigned int n_samples,
               unsigned int stripe_start,
               unsigned int stripe_stop,
               bool bypass_tips,
               unsigned int nthreads);

#endif
