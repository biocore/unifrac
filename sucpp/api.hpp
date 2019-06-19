#include "task_parameters.hpp"

#ifdef __cplusplus
#include <vector>
#define EXTERN extern "C"


#else
#include <stdbool.h>
#define EXTERN
#endif

#define PARTIAL_MAGIC "SSU-PARTIAL-01"

typedef enum compute_status {okay=0, tree_missing, table_missing, table_empty, unknown_method, table_and_tree_do_not_overlap} ComputeStatus;
typedef enum io_status {read_okay=0, write_okay, open_error, read_error, magic_incompatible, bad_header, unexpected_end} IOStatus;
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

void destroy_mat(mat_t** result);
void destroy_partial_mat(partial_mat_t** result);
void destroy_results_vec(r_vec** result);

/* Compute UniFrac
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
EXTERN IOStatus write_partial(const char* filename, partial_mat_t* result);

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

/* Merge partial results
 *
 * results <partial_mat_t**> an array of partial_mat_t*
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
