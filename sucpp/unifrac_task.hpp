#include <math.h>
#include <vector>
#include <stdint.h>

namespace su {
    struct task_parameters {
       uint32_t n_samples;          // number of samples
       unsigned int start;          // starting stripe
       unsigned int stop;           // stopping stripe
       unsigned int tid;            // thread ID
       double g_unifrac_alpha;      // generalized unifrac alpha
    };

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
