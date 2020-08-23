#include <algorithm> 
#include "unifrac_task.hpp"
#include <cstdlib>




template<class TFloat>
void su::UnifracUnnormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;

    const unsigned int step_size = this->step_size;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,dm_stripes_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx)*n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];

            // no need, if step_size<= UNIFRAC_BLOCK 
            if (k>=n_samples) continue; // past the limit

            const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

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

#ifdef _OPENACC
   // next iteration will use the alternative space
   std::swap(this->embedded_proportions,this->embedded_proportions_alt);
#endif
}

template<class TFloat>
void su::UnifracVawUnnormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
    const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;

    const unsigned int step_size = this->step_size;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];

            if (k>=n_samples) continue; // past the limit

            const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r*emb;

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

#ifdef _OPENACC
   // next iteration will use the alternative space
   std::swap(this->embedded_proportions,this->embedded_proportions_alt);
#endif
}

template<class TFloat>
void su::UnifracNormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;

    const unsigned int step_size = this->step_size;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
 	for(unsigned int ik = 0; ik < step_size ; ik++) {
	    const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

	    // no need if step_size<=UNIFRAC_BLOCK 
	    if (k>=n_samples) continue; // past the limit

            const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

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

#ifdef _OPENACC
   // next iteration will use the alternative space
   std::swap(this->embedded_proportions,this->embedded_proportions_alt);
#endif
}

template<class TFloat>
void su::UnifracVawNormalizedWeightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
    const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;

    const unsigned int step_size = this->step_size;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

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

#ifdef _OPENACC
   // next iteration will use the alternative space
   std::swap(this->embedded_proportions,this->embedded_proportions_alt);
#endif
}

template<class TFloat>
void su::UnifracGeneralizedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;

    const TFloat g_unifrac_alpha = this->task_p->g_unifrac_alpha;

    const unsigned int step_size = this->step_size;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

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

#ifdef _OPENACC
   // next iteration will use the alternative space
   std::swap(this->embedded_proportions,this->embedded_proportions_alt);
#endif
}

template<class TFloat>
void su::UnifracVawGeneralizedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;

    const TFloat g_unifrac_alpha = this->task_p->g_unifrac_alpha;

    // openacc only works well with local variables
    const TFloat * const __restrict__ embedded_proportions = this->embedded_proportions;
    const TFloat * const __restrict__ embedded_counts = this->embedded_counts;
    const TFloat * const __restrict__ sample_total_counts = this->sample_total_counts;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;

    const unsigned int step_size = this->step_size;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up
    // quick hack, to be finished

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

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

#ifdef _OPENACC
   // next iteration will use the alternative space
   std::swap(this->embedded_proportions,this->embedded_proportions_alt);
#endif
}

template<class TFloat>
void su::UnifracUnweightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;

    // openacc only works well with local variables
    const uint32_t * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;

    TFloat * const __restrict__ sums = this->sums;

    const unsigned int step_size = this->step_size;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    const unsigned int filled_embs_els = filled_embs/32;
    const unsigned int filled_embs_rem = filled_embs%32; 

    const unsigned int filled_embs_els_round = (filled_embs+31)/32;


    // pre-compute sums, since they are likely to be accessed many times
#pragma acc parallel loop collapse(2) gang present(lengths,sums) async
    for (unsigned int emb_el=0; emb_el<filled_embs_els; emb_el++) {
       for (unsigned int sub4=0; sub4<4; sub4++) {
          const unsigned int emb4 = emb_el*4+sub4;
          TFloat * __restrict__ psum = &(sums[emb4<<8]);
          const TFloat * __restrict__ pl   = &(lengths[emb4*8]);

#pragma acc loop vector
          // compute all the combinations for this block
          for (unsigned int b8_i=0; b8_i<0x100; b8_i++) {
             psum[b8_i] = (((b8_i >> 0) & 1) * pl[0]) + (((b8_i >> 1) & 1) * pl[1]) + 
                          (((b8_i >> 2) & 1) * pl[2]) + (((b8_i >> 3) & 1) * pl[3]) +
                          (((b8_i >> 4) & 1) * pl[4]) + (((b8_i >> 5) & 1) * pl[5]) +
                          (((b8_i >> 6) & 1) * pl[6]) + (((b8_i >> 7) & 1) * pl[7]);
          }
       }
    }
    if (filled_embs_rem>0) { // add also the overflow elements
       const unsigned int emb_el=filled_embs_els;
#pragma acc parallel loop gang present(lengths,sums) async
       for (unsigned int sub4=0; sub4<4; sub4++) {
          // we are summing we have enough buffer in sums
          const unsigned int emb4 = emb_el*4+sub4;
          TFloat * __restrict__ psum = &(sums[emb4<<8]);

#pragma acc loop vector
          // compute all the combinations for this block, set to 0 any past the limit
          for (unsigned int b8_i=0; b8_i<0x100; b8_i++) {
             TFloat val= 0;
             for (unsigned int li=(emb4*8); li<filled_embs; li++) {
               val += ((b8_i >>  (li-(emb4*8))) & 1) * lengths[li];
             }
             psum[b8_i] = val;
          }
        }
    }

    // point of thread
#ifdef _OPENACC
#pragma acc wait
#pragma acc parallel loop collapse(3) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,sums) async
#else
#pragma omp parallel for schedule(dynamic,1)
#endif
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#pragma acc loop seq
            for (unsigned int emb_el=0; emb_el<filled_embs_els_round; emb_el++) {
                const uint64_t offset = n_samples_r * emb_el;
                const TFloat * __restrict__ psum = &(sums[emb_el*0x400]);

                uint32_t u1 = embedded_proportions[offset + k];
                uint32_t v1 = embedded_proportions[offset + l1];
                uint32_t x1 = u1 ^ v1;
                uint32_t o1 = u1 | v1;

                // use the pre-computed sums

                my_stripe       += psum[              (x1 & 0xff)] + 
                                   psum[0x100+((x1 >>  8) & 0xff)] +
                                   psum[0x200+((x1 >> 16) & 0xff)] +
                                   psum[0x300+((x1 >> 24)       )];
                my_stripe_total += psum[              (o1 & 0xff)] +
                                   psum[0x100+((o1 >>  8) & 0xff)] +
                                   psum[0x200+((o1 >> 16) & 0xff)] +
                                   psum[0x300+((o1 >> 24)       )];
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }

#ifdef _OPENACC
   // next iteration will use the alternative space
   std::swap(this->embedded_proportions,this->embedded_proportions_alt);
#endif
}

template<class TFloat>
void su::UnifracVawUnweightedTask<TFloat>::_run(unsigned int filled_embs, const TFloat * __restrict__ lengths) {
    const unsigned int start_idx = this->task_p->start;
    const unsigned int stop_idx = this->task_p->stop;
    const unsigned int n_samples = this->task_p->n_samples;
    const uint64_t n_samples_r = this->dm_stripes.n_samples_r;

    // openacc only works well with local variables
    const uint32_t * const __restrict__ embedded_proportions = this->embedded_proportions;
    const TFloat  * const __restrict__ embedded_counts = this->embedded_counts;
    const TFloat  * const __restrict__ sample_total_counts = this->sample_total_counts;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;

    const unsigned int step_size = this->step_size;
    const unsigned int sample_steps = n_samples+(step_size-1)/step_size; // round up

    const unsigned int filled_embs_els = (filled_embs+31)/32; // round up

    // point of thread
#pragma acc parallel loop collapse(3) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
      for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
        for(unsigned int ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

#pragma acc loop seq
            for (unsigned int emb_el=0; emb_el<filled_embs_els; emb_el++) {
                const uint64_t offset_p = n_samples_r * emb_el;
                uint32_t u1 = embedded_proportions[offset_p + k];
                uint32_t v1 = embedded_proportions[offset_p + l1];
                uint32_t x1 = u1 ^ v1;
                uint32_t o1 = u1 | v1;

                // embedded_proporions is packed

#pragma acc loop seq
                for (unsigned int ei=0; ei<32; ei++) {
                   unsigned int emb=emb_el*32+ei;
                   if (emb<filled_embs) {
                     const uint64_t offset_c = n_samples_r * emb;
                     TFloat mi = embedded_counts[offset_c + k] + embedded_counts[offset_c + l1];
                     TFloat vaw = sqrt(mi * (m - mi));

                     // embedded_counts is not packed

                     if(vaw > 0) {
                       TFloat length = lengths[emb];
                       TFloat lv1 = length / vaw;

                       my_stripe +=       ((x1 >> ei) & 1)*lv1;
                       my_stripe_total += ((o1 >> ei) & 1)*lv1;
                     }
                   }
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }

#ifdef _OPENACC
   // next iteration will use the alternative space
   std::swap(this->embedded_proportions,this->embedded_proportions_alt);
#endif
}

