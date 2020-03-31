#include "unifrac_task.hpp"
#include <cstdlib>


void su::UnifracUnnormalizedWeightedTask::_run(unsigned int filled_embs, const double * restrict lengths) {
    const unsigned int start_idx = task_p->start;
    const unsigned int stop_idx = task_p->stop;
    const unsigned int n_samples = task_p->n_samples;
    const unsigned int trailing = n_samples - (n_samples % 4);

    // openacc only works well with local variables
    const double * const embedded_proportions = this->embedded_proportions;
    double * const dm_stripes_buf = this->dm_stripes.buf;

    // quick hack, to be finished
    const double length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples / 4; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            double *dm_stripe = dm_stripes_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];

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

#ifdef _OPENACC
    }

    if (trailing<n_samples) {
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf)
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++)
#endif
        for(unsigned int k = trailing; k < n_samples; k++) {
                unsigned int idx = (stripe-start_idx)*n_samples;
                double *dm_stripe = dm_stripes_buf+idx;
                //double *dm_stripe = dm_stripes[stripe];

                double u = embedded_proportions[k];
                double v = embedded_proportions[k + stripe + 1];
 
                dm_stripe[k] += fabs(u - v) * length;
        }
    }
}

void su::UnifracVawUnnormalizedWeightedTask::_run(unsigned int filled_embs, const double * restrict lengths) {
    const unsigned int start_idx = task_p->start;
    const unsigned int stop_idx = task_p->stop;
    const unsigned int n_samples = task_p->n_samples;

    // openacc only works well with local variables
    const double * const embedded_proportions = this->embedded_proportions;
    const double * const embedded_counts = this->embedded_counts;
    const double * const sample_total_counts = this->sample_total_counts;
    double * const dm_stripes_buf = this->dm_stripes.buf;

    // quick hack, to be finished
    const double length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            double *dm_stripe = dm_stripes_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];

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

void su::UnifracNormalizedWeightedTask::_run(unsigned int filled_embs, const double * restrict lengths) {
    const unsigned int start_idx = task_p->start;
    const unsigned int stop_idx = task_p->stop;
    const unsigned int n_samples = task_p->n_samples;

    // openacc only works well with local variables
    const double * const restrict embedded_proportions = this->embedded_proportions;
    double * const restrict dm_stripes_buf = this->dm_stripes.buf;
    double * const restrict dm_stripes_total_buf = this->dm_stripes_total.buf;


    const unsigned int step_size = 128;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
 	for(unsigned int ik = 0; ik < step_size ; ik++) {
	    uint64_t k = sk*step_size + ik;
            uint64_t idx = (stripe-start_idx);
            idx *= n_samples; // force 64 bit multiply
            double * restrict dm_stripe = dm_stripes_buf+idx;
            double * restrict dm_stripe_total = dm_stripes_total_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];
            //double *dm_stripe_total = dm_stripes_total[stripe];

            uint64_t l = k + stripe;

	    if (k>=n_samples) continue; // past the limit

            double my_stripe = dm_stripe[k];
            double my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                uint64_t offset = n_samples*2;
                offset *= emb; // force 64-bit multiply

                double u1 = embedded_proportions[offset + k];
                double v1 = embedded_proportions[offset + l + 1];
                double diff1 = u1 - v1;
                double sum1 = u1 + v1;
                double length = lengths[emb];

                my_stripe     += fabs(diff1) * length;
                my_stripe_total     += sum1 * length;
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }

}

void su::UnifracVawNormalizedWeightedTask::_run(unsigned int filled_embs, const double * restrict lengths) {
    const unsigned int start_idx = task_p->start;
    const unsigned int stop_idx = task_p->stop;
    const unsigned int n_samples = task_p->n_samples;

    // openacc only works well with local variables
    const double * const embedded_proportions = this->embedded_proportions;
    const double * const embedded_counts = this->embedded_counts;
    const double * const sample_total_counts = this->sample_total_counts;
    double * const dm_stripes_buf = this->dm_stripes.buf;
    double * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const double length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            double *dm_stripe = dm_stripes_buf+idx;
            double *dm_stripe_total = dm_stripes_total_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];
            //double *dm_stripe_total = dm_stripes_total[stripe];

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
                                   double sum_pow1 = pow(s, g_unifrac_alpha) * length; \
                                   dm_stripe[j] += sum_pow1 * (sub1 / s); \
                                   dm_stripe_total[j] += sum_pow1; \
                               }
void su::UnifracGeneralizedTask::_run(unsigned int filled_embs, const double * restrict lengths) {
    const unsigned int start_idx = task_p->start;
    const unsigned int stop_idx = task_p->stop;
    const unsigned int n_samples = task_p->n_samples;
    const unsigned int trailing = n_samples - (n_samples % 4);

    const double g_unifrac_alpha = task_p->g_unifrac_alpha;

    // openacc only works well with local variables
    const double * const embedded_proportions = this->embedded_proportions;
    double * const dm_stripes_buf = this->dm_stripes.buf;
    double * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const double length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples / 4; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            double *dm_stripe = dm_stripes_buf+idx;
            double *dm_stripe_total = dm_stripes_total_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];
            //double *dm_stripe_total = dm_stripes_total[stripe];

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
        

#ifdef _OPENACC
    }

    if (trailing<n_samples) {
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf)
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++)
#endif
        for(unsigned int k = trailing; k < n_samples; k++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            double *dm_stripe = dm_stripes_buf+idx;
            double *dm_stripe_total = dm_stripes_total_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];
            //double *dm_stripe_total = dm_stripes_total[stripe];

            double u = embedded_proportions[k];
            double v = embedded_proportions[k + stripe + 1];
            double s = u + v;
            GUNIFRAC(u, v, s, k)
        }
    }
}

void su::UnifracVawGeneralizedTask::_run(unsigned int filled_embs, const double * restrict lengths) {
    const unsigned int start_idx = task_p->start;
    const unsigned int stop_idx = task_p->stop;
    const unsigned int n_samples = task_p->n_samples;

    const double g_unifrac_alpha = task_p->g_unifrac_alpha;

    // openacc only works well with local variables
    const double * const embedded_proportions = this->embedded_proportions;
    const double * const embedded_counts = this->embedded_counts;
    const double * const sample_total_counts = this->sample_total_counts;
    double * const dm_stripes_buf = this->dm_stripes.buf;
    double * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const double length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            double *dm_stripe = dm_stripes_buf+idx;
            double *dm_stripe_total = dm_stripes_total_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];
            //double *dm_stripe_total = dm_stripes_total[stripe];

            double m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            double mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            double vaw = sqrt(mi * (m - mi));
            
            double u1 = embedded_proportions[j];
            double v1 = embedded_proportions[j + stripe + 1];
            
            if(vaw > 0.0) {
                double sum1 = (u1 + v1) / vaw;
                double sub1 = fabs(u1 - v1) / vaw;
                double sum_pow1 = pow(sum1, g_unifrac_alpha) * length;
                dm_stripe[j] += sum_pow1 * (sub1 / sum1);
                dm_stripe_total[j] += sum_pow1;
            }
        }
    }
}
void su::UnifracUnweightedTask::_run(unsigned int filled_embs, const double * restrict lengths) {
    const unsigned int start_idx = task_p->start;
    const unsigned int stop_idx = task_p->stop;
    const unsigned int n_samples = task_p->n_samples;
    const unsigned int trailing = n_samples - (n_samples % 4);

    // openacc only works well with local variables
    const double * const embedded_proportions = this->embedded_proportions;
    double * const dm_stripes_buf = this->dm_stripes.buf;
    double * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const double length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples / 4; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            double *dm_stripe = dm_stripes_buf+idx;
            double *dm_stripe_total = dm_stripes_total_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];
            //double *dm_stripe_total = dm_stripes_total[stripe];

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
       
#ifdef _OPENACC
    }

    if (trailing<n_samples) {
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf)
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++)
#endif
        for(unsigned int k = trailing; k < n_samples; k++) {
                unsigned int idx = (stripe-start_idx)*n_samples;
                double *dm_stripe = dm_stripes_buf+idx;
                double *dm_stripe_total = dm_stripes_total_buf+idx;
                //double *dm_stripe = dm_stripes[stripe];
                //double *dm_stripe_total = dm_stripes_total[stripe];

                int32_t u = embedded_proportions[k] > 0;
                int32_t v = embedded_proportions[k + stripe + 1] > 0;

                dm_stripe[k] += (u ^ v) * length;
                dm_stripe_total[k] += (u | v) * length;
        }
    }
}

void su::UnifracVawUnweightedTask::_run(unsigned int filled_embs, const double * restrict lengths) {
    const unsigned int start_idx = task_p->start;
    const unsigned int stop_idx = task_p->stop;
    const unsigned int n_samples = task_p->n_samples;

    // openacc only works well with local variables
    const double * const embedded_proportions = this->embedded_proportions;
    const double * const embedded_counts = this->embedded_counts;
    const double * const sample_total_counts = this->sample_total_counts;
    double * const dm_stripes_buf = this->dm_stripes.buf;
    double * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const double length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            double *dm_stripe = dm_stripes_buf+idx;
            double *dm_stripe_total = dm_stripes_total_buf+idx;
            //double *dm_stripe = dm_stripes[stripe];
            //double *dm_stripe_total = dm_stripes_total[stripe];

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

