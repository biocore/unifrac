#include "unifrac_task.hpp"
#include <cstdlib>


void su::_unnormalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                             std::vector<double*> &__restrict__ dm_stripes_total,
                                             double* __restrict__ embedded_proportions,
                                             double length,
                                             const su::task_parameters* task_p) {
    double *dm_stripe;
    for(unsigned int stripe=task_p->start; stripe < task_p->stop; stripe++) {
        dm_stripe = dm_stripes[stripe];

        /* intrinsics yield about a 2x reduction in runtime on llvm. they
         * were not effective on linux gcc 4.9.1 or 4.9.2. it is unclear
         * if they would be effective on other versions of gcc.
         *
         * one reason they help is that these for loops are not easily
         * autovectorizable. using the intrinsics effectively gets around
         * this. ...although, it also appears that loop unrolling works.
         *
         * it may make sense to revisit the inclusion of intriniscs, however
         * support must be tested at compile time, so it's rather annoying
         * at the moment. basically, we can't assume the presence of avx2.
         */
        for(unsigned int j = 0; j < task_p->n_samples / 4; j++) {
            int k = j * 4;
            double u1 = embedded_proportions[k];
            double u2 = embedded_proportions[k + 1];
            double u3 = embedded_proportions[k + 2];
            double u4 = embedded_proportions[k + 3];

            double v1 = embedded_proportions[k + stripe + 1];
            double v2 = embedded_proportions[k + stripe + 2];
            double v3 = embedded_proportions[k + stripe + 3];
            double v4 = embedded_proportions[k + stripe + 4];
            
            dm_stripe[k] += fabs(u1 - v1) * length;
            dm_stripe[k + 1] += fabs(u2 - v2) * length;
            dm_stripe[k + 2] += fabs(u3 - v3) * length;
            dm_stripe[k + 3] += fabs(u4 - v4) * length;
        }

        if((task_p->n_samples % 4) != 0) {
            for(unsigned int k = task_p->n_samples - (task_p->n_samples % 4); k < task_p->n_samples; k++) {
                double u = embedded_proportions[k];
                double v = embedded_proportions[k + stripe + 1];
 
                dm_stripe[k] += fabs(u - v) * length;
            }
        }
    }
}

void su::_vaw_unnormalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                                 std::vector<double*> &__restrict__ dm_stripes_total,
                                                 double* __restrict__ embedded_proportions,
                                                 double* __restrict__ embedded_counts,
                                                 double* __restrict__ sample_total_counts,
                                                 double length,
                                                 const su::task_parameters* task_p) {
    double *dm_stripe;
    for(unsigned int stripe=task_p->start; stripe < task_p->stop; stripe++) {
        dm_stripe = dm_stripes[stripe];

        for(unsigned int j = 0; j < task_p->n_samples; j++) {
            double u = embedded_proportions[j];
            double v = embedded_proportions[j + stripe + 1];

            double m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            double mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            double vaw = sqrt(mi * (m - mi));

            if(vaw > 0)
                dm_stripe[j] += (fabs(u - v) * length) / vaw;
        }
    }
}
void su::_normalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                           std::vector<double*> &__restrict__ dm_stripes_total,
                                           double* __restrict__ embedded_proportions, 
                                           double length, 
                                           const su::task_parameters* task_p) {
    double *dm_stripe;
    double *dm_stripe_total;
    unsigned int trailing = task_p->n_samples - (task_p->n_samples % 4);

    // point of thread
    for(unsigned int stripe = task_p->start; stripe < task_p->stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < task_p->n_samples / 4; j++) {
            int k = j * 4;
            int l = k + stripe;

            double u1 = embedded_proportions[k];
            double u2 = embedded_proportions[k + 1];
            double u3 = embedded_proportions[k + 2];
            double u4 = embedded_proportions[k + 3];
         
            double v1 = embedded_proportions[l + 1];
            double v2 = embedded_proportions[l + 2];
            double v3 = embedded_proportions[l + 3];
            double v4 = embedded_proportions[l + 4];
            
            double diff1 = u1 - v1;   
            double diff2 = u2 - v2;   
            double diff3 = u3 - v3;   
            double diff4 = u4 - v4;

            double sum1 = u1 + v1;   
            double sum2 = u2 + v2;   
            double sum3 = u3 + v3;   
            double sum4 = u4 + v4;

            dm_stripe[k]     += fabs(diff1) * length;
            dm_stripe[k + 1] += fabs(diff2) * length;
            dm_stripe[k + 2] += fabs(diff3) * length;
            dm_stripe[k + 3] += fabs(diff4) * length;

            dm_stripe_total[k]     += sum1 * length;
            dm_stripe_total[k + 1] += sum2 * length;
            dm_stripe_total[k + 2] += sum3 * length;
            dm_stripe_total[k + 3] += sum4 * length;
        }

        for(unsigned int k = trailing; k < task_p->n_samples; k++) {
            double u = embedded_proportions[k];
            double v = embedded_proportions[k + stripe + 1];
            double diff = u - v;   
            double sum = u + v;   
               
            dm_stripe[k] += fabs(diff) * length;
            dm_stripe_total[k] += sum * length;
        }
    }
}

void su::_vaw_normalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                               std::vector<double*> &__restrict__ dm_stripes_total,
                                               double* __restrict__ embedded_proportions, 
                                               double* __restrict__ embedded_counts, 
                                               double* __restrict__ sample_total_counts,
                                               double length, 
                                               const su::task_parameters* task_p) {
    double *dm_stripe;
    double *dm_stripe_total;

    // point of thread
    for(unsigned int stripe = task_p->start; stripe < task_p->stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < task_p->n_samples; j++) {
            double u = embedded_proportions[j];
            double v = embedded_proportions[j + stripe + 1];
            
            double m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            double mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            double vaw = sqrt(mi * (m - mi));

            if(vaw > 0) {   
                dm_stripe[j] += (fabs(u - v) * length) / vaw;
                dm_stripe_total[j] += ((u + v) * length) / vaw;
            }
        }
    }
}

#define GUNIFRAC(u, v, s, j)   if(s != 0.0) { \
                                   double sub1 = fabs(u - v); \
                                   double sum_pow1 = pow(s, task_p->g_unifrac_alpha) * length; \
                                   dm_stripe[j] += sum_pow1 * (sub1 / s); \
                                   dm_stripe_total[j] += sum_pow1; \
                               }
void su::_generalized_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                   std::vector<double*> &__restrict__ dm_stripes_total,
                                   double* __restrict__ embedded_proportions, 
                                   double length, 
                                   const su::task_parameters* task_p) {
    double *dm_stripe;
    double *dm_stripe_total;
    unsigned int trailing = task_p->n_samples - (task_p->n_samples % 4);

    // point of thread
    for(unsigned int stripe = task_p->start; stripe < task_p->stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < task_p->n_samples / 4; j++) {
            int k = j * 4;
            int l = k + stripe;

            double u1 = embedded_proportions[k];
            double u2 = embedded_proportions[k + 1];
            double u3 = embedded_proportions[k + 2];
            double u4 = embedded_proportions[k + 3];
         
            double v1 = embedded_proportions[l + 1];
            double v2 = embedded_proportions[l + 2];
            double v3 = embedded_proportions[l + 3];
            double v4 = embedded_proportions[l + 4];
            
            double sum1 = u1 + v1;
            double sum2 = u2 + v2;
            double sum3 = u3 + v3;
            double sum4 = u4 + v4;

            GUNIFRAC(u1, v1, sum1, k)
            GUNIFRAC(u2, v2, sum2, k + 1)
            GUNIFRAC(u3, v3, sum3, k + 2)
            GUNIFRAC(u4, v4, sum4, k + 3)
        }
        
        for(unsigned int k = trailing; k < task_p->n_samples; k++) {
            double u = embedded_proportions[k];
            double v = embedded_proportions[k + stripe + 1];
            double s = u + v;
            GUNIFRAC(u, v, s, k)
        }
    }
}

void su::_vaw_generalized_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                       std::vector<double*> &__restrict__ dm_stripes_total,
                                       double* __restrict__ embedded_proportions, 
                                       double* __restrict__ embedded_counts, 
                                       double* __restrict__ sample_total_counts,
                                       double length, 
                                       const su::task_parameters* task_p) {
    double *dm_stripe;
    double *dm_stripe_total;

    // point of thread
    for(unsigned int stripe = task_p->start; stripe < task_p->stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < task_p->n_samples; j++) {
            double m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            double mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            double vaw = sqrt(mi * (m - mi));
            
            double u1 = embedded_proportions[j];
            double v1 = embedded_proportions[j + stripe + 1];
            
            if(vaw > 0.0) {
                double sum1 = (u1 + v1) / vaw;
                double sub1 = fabs(u1 - v1) / vaw;
                double sum_pow1 = pow(sum1, task_p->g_unifrac_alpha) * length;
                dm_stripe[j] += sum_pow1 * (sub1 / sum1);
                dm_stripe_total[j] += sum_pow1;
            }
        }
    }
}
void su::_unweighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                  std::vector<double*> &__restrict__ dm_stripes_total,
                                  double* __restrict__ embedded_proportions, 
                                  double length,  
                                  const su::task_parameters* task_p) {
    double *dm_stripe;
    double *dm_stripe_total;
    
    for(unsigned int stripe = task_p->start; stripe < task_p->stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < task_p->n_samples / 4; j++) {
            int k = j * 4;
            int32_t u1 = embedded_proportions[k] > 0;
            int32_t u2 = embedded_proportions[k + 1] > 0;
            int32_t u3 = embedded_proportions[k + 2] > 0;
            int32_t u4 = embedded_proportions[k + 3] > 0;
 
            int32_t v1 = embedded_proportions[k + stripe + 1] > 0;
            int32_t v2 = embedded_proportions[k + stripe + 2] > 0;
            int32_t v3 = embedded_proportions[k + stripe + 3] > 0;
            int32_t v4 = embedded_proportions[k + stripe + 4] > 0;

            dm_stripe[k] += (u1 ^ v1) * length;
            dm_stripe[k + 1] += (u2 ^ v2) * length;
            dm_stripe[k + 2] += (u3 ^ v3) * length;
            dm_stripe[k + 3] += (u4 ^ v4) * length;
 
            dm_stripe_total[k] += (u1 | v1) * length;
            dm_stripe_total[k + 1] += (u2 | v2) * length;
            dm_stripe_total[k + 2] += (u3 | v3) * length;
            dm_stripe_total[k + 3] += (u4 | v4) * length;
        }
        
        if((task_p->n_samples % 4) != 0) {
            for(unsigned int k = task_p->n_samples - (task_p->n_samples % 4); k < task_p->n_samples; k++) {
                int32_t u = embedded_proportions[k] > 0;
                int32_t v = embedded_proportions[k + stripe + 1] > 0;

                dm_stripe[k] += (u ^ v) * length;
                dm_stripe_total[k] += (u | v) * length;
            }
        }
    }
}

void su::_vaw_unweighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                      std::vector<double*> &__restrict__ dm_stripes_total,
                                      double* __restrict__ embedded_proportions, 
                                      double* __restrict__ embedded_counts, 
                                      double* __restrict__ sample_total_counts,
                                      double length,  
                                      const su::task_parameters* task_p) {
    double *dm_stripe;
    double *dm_stripe_total;
    
    for(unsigned int stripe = task_p->start; stripe < task_p->stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < task_p->n_samples; j++) {
            int32_t u = embedded_proportions[j] > 0;
            int32_t v = embedded_proportions[j + stripe + 1] > 0;

            double m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            double mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            double vaw = sqrt(mi * (m - mi));
            
            if(vaw > 0) {
                dm_stripe[j] += ((u ^ v) * length) / vaw;
                dm_stripe_total[j] += ((u | v) * length) / vaw;
            }
        }
    }
}

