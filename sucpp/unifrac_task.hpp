#include <math.h>
#include <vector>
#include <stdint.h>

namespace su {
    /* task specific compute parameters
     *
     * n_samples <int> the number of samples being processed
     * start <uint> the first stride to process
     * stop <uint> the last stride to process
     * tid <uint> the thread identifier
     * g_unifrac_alpha <double> an alpha value for generalized unifrac
     */
    struct task_parameters {
       uint32_t n_samples;          // number of samples
       unsigned int start;          // starting stripe
       unsigned int stop;           // stopping stripe
       unsigned int tid;            // thread ID

       // task specific arguments below
       double g_unifrac_alpha;      // generalized unifrac alpha
    };

    /* void su::unifrac tasks
     *
     * all methods utilize the same function signature. that signature is as follows:
     *
     * dm_stripes vector<double> the stripes of the distance matrix being accumulated 
     *      into for unique branch length
     * dm_stripes vector<double> the stripes of the distance matrix being accumulated 
     *      into for total branch length (e.g., to normalize unweighted unifrac)
     * embedded_proportions <double*> the proportions vector for a sample, or rather
     *      the counts vector normalized to 1. this vector is embedded as it is 
     *      duplicated: if A, B and C are proportions for features A, B, and C, the
     *      vector will look like [A B C A B C].
     * length <double> the branch length of the current node to its parent.
     * tasp_p <task_parameters*> task specific parameters.
    void _unnormalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                             std::vector<double*> &__restrict__ dm_stripes_total,
                                             double* __restrict__ embedded_proportions,
                                             double length,
                                             const su::task_parameters* task_p);
    void _normalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                           std::vector<double*> &__restrict__ dm_stripes_total,
                                           double* __restrict__ embedded_proportions,
                                           double length,
                                           const su::task_parameters* task_p);
    void _unweighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                  std::vector<double*> &__restrict__ dm_stripes_total,
                                  double* __restrict__ embedded_proportions,
                                  double length,
                                  const su::task_parameters* task_p);
    void _generalized_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                   std::vector<double*> &__restrict__ dm_stripes_total,
                                   double* __restrict__ embedded_proportions,
                                   double length,
                                   const su::task_parameters* task_p);
    
    /* void su::unifrac_vaw tasks
     *
     * all methods utilize the same function signature. that signature is as follows:
     *
     * dm_stripes vector<double> the stripes of the distance matrix being accumulated 
     *      into for unique branch length
     * dm_stripes vector<double> the stripes of the distance matrix being accumulated 
     *      into for total branch length (e.g., to normalize unweighted unifrac)
     * embedded_proportions <double*> the proportions vector for a sample, or rather
     *      the counts vector normalized to 1. this vector is embedded as it is 
     *      duplicated: if A, B and C are proportions for features A, B, and C, the
     *      vector will look like [A B C A B C].
     * embedded_counts <double*> the counts vector embedded in the same way and order as
     *      embedded_proportions. the values of this array are unnormalized feature 
     *      counts for the subtree.
     * sample_total_counts <double*> the total unnormalized feature counts for all samples
     *      embedded in the same way and order as embedded_proportions.
     * length <double> the branch length of the current node to its parent.
     * tasp_p <task_parameters*> task specific parameters.
    void _vaw_unnormalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                                 std::vector<double*> &__restrict__ dm_stripes_total,
                                                 double* __restrict__ embedded_proportions,
                                                 double* __restrict__ embedded_counts,
                                                 double* __restrict__ sample_total_counts,
                                                 double length,
                                                 const su::task_parameters* task_p);
    void _vaw_normalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                               std::vector<double*> &__restrict__ dm_stripes_total,
                                               double* __restrict__ embedded_proportions,
                                               double* __restrict__ embedded_counts,
                                               double* __restrict__ sample_total_counts,
                                               double length,
                                               const su::task_parameters* task_p);
    void _vaw_unweighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                      std::vector<double*> &__restrict__ dm_stripes_total,
                                      double* __restrict__ embedded_proportions,
                                      double* __restrict__ embedded_counts,
                                      double* __restrict__ sample_total_counts,
                                      double length,
                                      const su::task_parameters* task_p);
    void _vaw_generalized_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                       std::vector<double*> &__restrict__ dm_stripes_total,
                                       double* __restrict__ embedded_proportions,
                                       double* __restrict__ embedded_counts,
                                       double* __restrict__ sample_total_counts,
                                       double length,
                                       const su::task_parameters* task_p);
}
