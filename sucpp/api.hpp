#include "task_parameters.hpp"

#ifdef __cplusplus
#include <vector>
#define EXTERN extern "C"


#else
#include <stdbool.h>
#define EXTERN 
#endif

typedef enum compute_status {okay, tree_missing, table_missing, unknown_method} ComputeStatus;

/* a result matrix
 *
 * n_samples <uint> the number of samples.
 * cf_size <uint> the size of the condensed form.
 * is_square <bool> if true, indicates condensed_form represents a square
 *      matrix, and only the upper triangle is contained. if false, 
 *      condensed_form represents a dissimilarity matrix in which case 
 *      condensed_form stores all non-diagonal values in row order.
 * condensed_form <double*> the matrix values of length cf_size.
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct mat {
    unsigned int n_samples;  
    unsigned int cf_size;   
    bool is_square;          
    double* condensed_form;    
    char** sample_ids;
} mat_t;

/* a partial result containing stripe data
 *
 * n_samples <uint> the number of samples.
 * sample_ids <char**> the sample IDs of length n_samples.
 * stripes <double**> the stripe data of dimension (stripe_stop - stripe_start, n_samples)
 * stripe_start <uint> the logical starting stripe in the final matrix. 
 * stripe_Stop <uint> the logical stopping stripe in the final matrix.
 */
typedef struct partial_mat {
    uint32_t n_samples;  
    char** sample_ids;
    double** stripes;
    uint32_t stripe_start;
    uint32_t stripe_stop;
} partial_mat_t;

/* Compute UniFrac
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * threads <uint> the number of threads to use.
 * result <mat*> the resulting distance matrix in condensed form
 *
 * one_off returns the following error codes:
 *
 * OK             : no problems encountered
 * TABLE_MISSING  : the filename for the table does not exist
 * TREE_MISSING   : the filename for the tree does not exist
 * UNKNOWN_METHOD : the requested method is unknown.
 */
EXTERN ComputeStatus one_off(const char* biom_filename, const char* tree_filename, 
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             unsigned int threads, mat_t** result);

void destroy_mat(mat_t** result);
void destroy_partial_mat(partial_mat_t** result);

/* Compute a subset of a UniFrac distance matrix
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * threads <uint> the number of threads to use.
 * stripe_start <uint> the starting stripe to compute
 * stripe_stop <uint> the last stripe to compute
 * dm_stripes <vector of double*> the unique branch length stripes. This is expected to be
 *      uninitialized, and is an output parameter.
 * dm_stripes_total <vector of double*> the total branch length stripes. This is expected to be
 *      uninitialized, and is an output parameter.
 *
 * partial returns the following error codes:
 *
 * OK             : no problems encountered
 * TABLE_MISSING  : the filename for the table does not exist
 * TREE_MISSING   : the filename for the tree does not exist
 * UNKNOWN_METHOD : the requested method is unknown.
 */

compute_status partial(const char* biom_filename, const char* tree_filename, 
                       const char* unifrac_method, bool variance_adjust, double alpha,
                       unsigned int threads, unsigned int stripe_start, unsigned int stripe_stop,
                       partial_mat_t** result);

#ifdef __cplusplus
// TODO: only needed for testing, should be encased in a macro
void set_tasks(std::vector<su::task_parameters> &tasks,
               double alpha,
               unsigned int n_samples,
               unsigned int stripe_start, 
               unsigned int stripe_stop, 
               unsigned int nthreads);
#endif
