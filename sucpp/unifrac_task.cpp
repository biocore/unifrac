#include "unifrac_task.hpp"
#include <cstdlib>


template<class TFloat>
void su::UnifracUnnormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const unsigned int trailing = n_samples - (n_samples % 4);

    // openacc only works well with local variables
    const TFloat * const embedded_proportions = this->embedded_proportions;
    TFloat * const dm_stripes_buf = this->dm_stripes.buf;

    // quick hack, to be finished
    const TFloat length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples / 4; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            TFloat *dm_stripe = dm_stripes_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];

            int k = j * 4;
            TFloat u1 = embedded_proportions[k];
            TFloat u2 = embedded_proportions[k + 1];
            TFloat u3 = embedded_proportions[k + 2];
            TFloat u4 = embedded_proportions[k + 3];

            TFloat v1 = embedded_proportions[k + stripe + 1];
            TFloat v2 = embedded_proportions[k + stripe + 2];
            TFloat v3 = embedded_proportions[k + stripe + 3];
            TFloat v4 = embedded_proportions[k + stripe + 4];
            
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
                TFloat *dm_stripe = dm_stripes_buf+idx;
                //TFloat *dm_stripe = dm_stripes[stripe];

                TFloat u = embedded_proportions[k];
                TFloat v = embedded_proportions[k + stripe + 1];
 
                dm_stripe[k] += fabs(u - v) * length;
        }
    }
}

template<class TFloat>
void su::UnifracVawUnnormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    // openacc only works well with local variables
    const TFloat * const embedded_proportions = this->embedded_proportions;
    const TFloat * const embedded_counts = this->embedded_counts;
    const TFloat * const sample_total_counts = this->sample_total_counts;
    TFloat * const dm_stripes_buf = this->dm_stripes.buf;

    // quick hack, to be finished
    const TFloat length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            TFloat *dm_stripe = dm_stripes_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];

            TFloat u = embedded_proportions[j];
            TFloat v = embedded_proportions[j + stripe + 1];

            TFloat m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            TFloat mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            TFloat vaw = sqrt(mi * (m - mi));

            if(vaw > 0)
                dm_stripe[j] += (fabs(u - v) * length) / vaw;
        }
    }
}

template<class TFloat>
void su::UnifracNormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;


#ifdef _OPENACC
    // The parallel nature of GPUs needs a largish step
    const unsigned int step_size = 16;
#else
    // The serial nature of CPU cores prefers a small step
    const unsigned int step_size = 4;
#endif
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
 	for(unsigned int ik = 0; ik < step_size ; ik++) {
	    uint64_t k = sk*step_size + ik;
            uint64_t idx = (stripe-start_idx);
            idx *= n_samples; // force 64 bit multiply
            TFloat * __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

	    if (k>=n_samples) continue; // past the limit

            unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                uint64_t offset = n_samples;
                offset *= emb; // force 64-bit multiply

                TFloat u1 = embedded_proportions[offset + k];
                TFloat v1 = embedded_proportions[offset + l1];
                TFloat diff1 = u1 - v1;
                TFloat sum1 = u1 + v1;
                TFloat length = lengths[emb];

                my_stripe     += fabs(diff1) * length;
                my_stripe_total     += sum1 * length;
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }

}

template<class TFloat>
void su::UnifracVawNormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    // openacc only works well with local variables
    const TFloat * const embedded_proportions = this->embedded_proportions;
    const TFloat * const embedded_counts = this->embedded_counts;
    const TFloat * const sample_total_counts = this->sample_total_counts;
    TFloat * const dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const TFloat length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            TFloat *dm_stripe = dm_stripes_buf+idx;
            TFloat *dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            TFloat u = embedded_proportions[j];
            TFloat v = embedded_proportions[j + stripe + 1];
            
            TFloat m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            TFloat mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            TFloat vaw = sqrt(mi * (m - mi));

            if(vaw > 0) {   
                dm_stripe[j] += (fabs(u - v) * length) / vaw;
                dm_stripe_total[j] += ((u + v) * length) / vaw;
            }
        }
    }
}

#define GUNIFRAC(u, v, s, j)   if(s != 0.0) { \
                                   TFloat sub1 = fabs(u - v); \
                                   TFloat sum_pow1 = pow(s, g_unifrac_alpha) * length; \
                                   dm_stripe[j] += sum_pow1 * (sub1 / s); \
                                   dm_stripe_total[j] += sum_pow1; \
                               }
template<class TFloat>
void su::UnifracGeneralizedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const unsigned int trailing = n_samples - (n_samples % 4);

    const TFloat g_unifrac_alpha = this->task_p->g_unifrac_alpha;

    // openacc only works well with local variables
    const TFloat * const embedded_proportions = this->embedded_proportions;
    TFloat * const dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const TFloat length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples / 4; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            TFloat *dm_stripe = dm_stripes_buf+idx;
            TFloat *dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            int k = j * 4;
            int l = k + stripe;

            TFloat u1 = embedded_proportions[k];
            TFloat u2 = embedded_proportions[k + 1];
            TFloat u3 = embedded_proportions[k + 2];
            TFloat u4 = embedded_proportions[k + 3];
         
            TFloat v1 = embedded_proportions[l + 1];
            TFloat v2 = embedded_proportions[l + 2];
            TFloat v3 = embedded_proportions[l + 3];
            TFloat v4 = embedded_proportions[l + 4];
            
            TFloat sum1 = u1 + v1;
            TFloat sum2 = u2 + v2;
            TFloat sum3 = u3 + v3;
            TFloat sum4 = u4 + v4;

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
            TFloat *dm_stripe = dm_stripes_buf+idx;
            TFloat *dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            TFloat u = embedded_proportions[k];
            TFloat v = embedded_proportions[k + stripe + 1];
            TFloat s = u + v;
            GUNIFRAC(u, v, s, k)
        }
    }
}

template<class TFloat>
void su::UnifracVawGeneralizedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    const TFloat g_unifrac_alpha = this->task_p->g_unifrac_alpha;

    // openacc only works well with local variables
    const TFloat * const embedded_proportions = this->embedded_proportions;
    const TFloat * const embedded_counts = this->embedded_counts;
    const TFloat * const sample_total_counts = this->sample_total_counts;
    TFloat * const dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const TFloat length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            TFloat *dm_stripe = dm_stripes_buf+idx;
            TFloat *dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            TFloat m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            TFloat mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            TFloat vaw = sqrt(mi * (m - mi));
            
            TFloat u1 = embedded_proportions[j];
            TFloat v1 = embedded_proportions[j + stripe + 1];
            
            if(vaw > 0.0) {
                TFloat sum1 = (u1 + v1) / vaw;
                TFloat sub1 = fabs(u1 - v1) / vaw;
                TFloat sum_pow1 = pow(sum1, g_unifrac_alpha) * length;
                dm_stripe[j] += sum_pow1 * (sub1 / sum1);
                dm_stripe_total[j] += sum_pow1;
            }
        }
    }
}
template<class TFloat>
void su::UnifracUnweightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const unsigned int trailing = n_samples - (n_samples % 4);

    // openacc only works well with local variables
    const TFloat * const embedded_proportions = this->embedded_proportions;
    TFloat * const dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const TFloat length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples / 4; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            TFloat *dm_stripe = dm_stripes_buf+idx;
            TFloat *dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

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
                TFloat *dm_stripe = dm_stripes_buf+idx;
                TFloat *dm_stripe_total = dm_stripes_total_buf+idx;
                //TFloat *dm_stripe = dm_stripes[stripe];
                //TFloat *dm_stripe_total = dm_stripes_total[stripe];

                int32_t u = embedded_proportions[k] > 0;
                int32_t v = embedded_proportions[k + stripe + 1] > 0;

                dm_stripe[k] += (u ^ v) * length;
                dm_stripe_total[k] += (u | v) * length;
        }
    }
}

template<class TFloat>
void su::UnifracVawUnweightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    // openacc only works well with local variables
    const TFloat * const embedded_proportions = this->embedded_proportions;
    const TFloat * const embedded_counts = this->embedded_counts;
    const TFloat * const sample_total_counts = this->sample_total_counts;
    TFloat * const dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const dm_stripes_total_buf = this->dm_stripes_total.buf;

    // quick hack, to be finished
    const TFloat length = lengths[0];

    // point of thread
#pragma acc parallel loop collapse(2) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf)
    for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int j = 0; j < n_samples; j++) {
            unsigned int idx = (stripe-start_idx)*n_samples;
            TFloat *dm_stripe = dm_stripes_buf+idx;
            TFloat *dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            int32_t u = embedded_proportions[j] > 0;
            int32_t v = embedded_proportions[j + stripe + 1] > 0;

            TFloat m = sample_total_counts[j] + sample_total_counts[j + stripe + 1];
            TFloat mi = embedded_counts[j] + embedded_counts[j + stripe + 1];
            TFloat vaw = sqrt(mi * (m - mi));
            
            if(vaw > 0) {
                dm_stripe[j] += ((u ^ v) * length) / vaw;
                dm_stripe_total[j] += ((u | v) * length) / vaw;
            }
        }
    }
}

