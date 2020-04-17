#include "unifrac_task.hpp"
#include <cstdlib>


template<class TFloat>
void su::UnifracUnnormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;

#ifdef _OPENACC
    // The parallel nature of GPUs needs a largish step
    const unsigned int step_size = 16;
#else
    // The serial nature of CPU cores prefers a small step
    const unsigned int step_size = 4;
#endif
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,dm_stripes_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            uint64_t k = sk*step_size + ik;
            uint64_t idx = (stripe-start_idx);
            idx *= n_samples; // force 64 bit multiply
            TFloat * __restrict__ dm_stripe = dm_stripes_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];

            if (k>=n_samples) continue; // past the limit

            unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                uint64_t offset = n_samples;
                offset *= emb; // force 64-bit multiply

                TFloat u1 = embedded_proportions[offset + k];
                TFloat v1 = embedded_proportions[offset + l1];
                TFloat diff1 = u1 - v1;
                TFloat length = lengths[emb];

                my_stripe     += fabs(diff1) * length;
            }

            dm_stripe[k]     = my_stripe;
        }

      }
    }
}

template<class TFloat>
void su::UnifracVawUnnormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
    const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;

#ifdef _OPENACC
    // The parallel nature of GPUs needs a largish step
    const unsigned int step_size = 16;
#else
    // The serial nature of CPU cores prefers a small step
    const unsigned int step_size = 4;
#endif
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            uint64_t k = sk*step_size + ik;
            uint64_t idx = (stripe-start_idx);
            idx *= n_samples; // force 64 bit multiply
            TFloat * __restrict__ dm_stripe = dm_stripes_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];

            if (k>=n_samples) continue; // past the limit

            unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                uint64_t offset = n_samples;
                offset *= emb; // force 64-bit multiply

                TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
                TFloat vaw = sqrt(mi * (m - mi));

                if(vaw > 0) {
                  TFloat u1 = embedded_proportions[offset + k];
                  TFloat v1 = embedded_proportions[offset + l1];
                  TFloat diff1 = fabs(u1 - v1);
                  TFloat length = lengths[emb];

                   my_stripe += (diff1 * length) / vaw;
                }
            }

            dm_stripe[k]     = my_stripe;
        }

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
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
    const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
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
#pragma acc parallel loop collapse(3) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
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

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                uint64_t offset = n_samples;
                offset *= emb; // force 64-bit multiply

                TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
                TFloat vaw = sqrt(mi * (m - mi));

                if(vaw > 0) {
                  TFloat u1 = embedded_proportions[offset + k];
                  TFloat v1 = embedded_proportions[offset + l1];
                  TFloat diff1 = fabs(u1 - v1);
                  TFloat length = lengths[emb];

                  my_stripe += (diff1 * length) / vaw;
                  my_stripe_total += ((u1 + v1) * length) / vaw;
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }
}

template<class TFloat>
void su::UnifracGeneralizedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;

    const TFloat g_unifrac_alpha = this->task_p->g_unifrac_alpha;

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
                TFloat sum1 = u1 + v1;

                if(sum1 != 0.0) { 
                   TFloat length = lengths[emb];
                   TFloat diff1 = fabs(u1 - v1);
                   TFloat sum_pow1 = pow(sum1, g_unifrac_alpha) * length; 

                   my_stripe += sum_pow1 * (diff1 / sum1); 
                   my_stripe_total += sum_pow1; 
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

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
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
    const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
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
    // quick hack, to be finished

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
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

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                uint64_t offset = n_samples;
                offset *= emb; // force 64-bit multiply

                TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
                TFloat vaw = sqrt(mi * (m - mi));

                if(vaw > 0) {
                  TFloat u1 = embedded_proportions[offset + k];
                  TFloat v1 = embedded_proportions[offset + l1];
                  TFloat length = lengths[emb];

                  TFloat sum1 = (u1 + v1) / vaw;
                  TFloat sub1 = fabs(u1 - v1) / vaw;
                  TFloat sum_pow1 = pow(sum1, g_unifrac_alpha) * length;

                  my_stripe += sum_pow1 * (sub1 / sum1);
                  my_stripe_total += sum_pow1;
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }
      }
    }
}

template<class TFloat>
void su::UnifracUnweightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
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

                int32_t u1 = embedded_proportions[offset + k]  > 0;
                int32_t v1 = embedded_proportions[offset + l1] > 0;
                int32_t x1 = u1 ^ v1;
                int32_t o1 = u1 | v1;
                TFloat length = lengths[emb];

                my_stripe     += x1 * length;
                my_stripe_total     += o1 * length;
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }
}

template<class TFloat>
void su::UnifracVawUnweightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
    const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
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
#pragma acc parallel loop collapse(3) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
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

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                uint64_t offset = n_samples;
                offset *= emb; // force 64-bit multiply

                TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
                TFloat vaw = sqrt(mi * (m - mi));

                if(vaw > 0) {
                  int32_t u = embedded_proportions[offset + k] > 0;
                  int32_t v = embedded_proportions[offset + l1] > 0;
                  TFloat length = lengths[emb];

                  my_stripe += ((u ^ v) * length) / vaw;
                  my_stripe_total += ((u | v) * length) / vaw;
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }
}

