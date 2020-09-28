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

    bool * const __restrict__ zcheck = this->zcheck;
    TFloat * const __restrict__ sums = this->sums;

    const unsigned int step_size = su::UnifracUnnormalizedWeightedTask<TFloat>::step_size;
    const unsigned int sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // check for zero values and pre-compute single column sums
#ifdef _OPENACC
#pragma acc parallel loop present(embedded_proportions,lengths,zcheck,sums)
#endif
    for(unsigned int k=0; k<n_samples; k++) {
            bool all_zeros=true;
            TFloat my_sum = 0.0;

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat u1 = embedded_proportions[offset + k];
                my_sum += u1*lengths[emb];
                all_zeros = all_zeros && (u1==0.0);
            }

            sums[k]     = my_sum;
            zcheck[k] = all_zeros;
    }


    // now do the real compute
#ifdef _OPENACC
    const unsigned int acc_vector_size = su::UnifracUnnormalizedWeightedTask<TFloat>::acc_vector_size;
#pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,dm_stripes_buf,lengths,zcheck,sums) async
#endif
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
     for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
      for(unsigned int ik = 0; ik < step_size ; ik++) {
       const unsigned int k = sk*step_size + ik;

       if (k>=n_samples) continue; // past the limit

       const bool zcheck_k = zcheck[sk]; // due to loop collapse in ACC, must load in here

       const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

       const bool allzero_k = zcheck[k];
       const bool allzero_l1 = zcheck[l1];

       if (allzero_k && allzero_l1) {
         // nothing to do, would have to add 0
       } else {
          TFloat my_stripe;

          if (allzero_k || allzero_l1) {
            // one side has all zeros
            // we can use the distributed property, and use the pre-computed values

            const unsigned int ridx = (allzero_k) ? l1 : // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
                                                    k;   // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0

            // keep reads in the same place to maximize GPU warp performance
            my_stripe = sums[ridx];

          } else {
            // both sides non zero, use the explicit but slow approach
            my_stripe = 0.0;

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat u1 = embedded_proportions[offset + k];
                TFloat v1 = embedded_proportions[offset + l1];
                TFloat diff1 = u1 - v1;
                TFloat length = lengths[emb];

                my_stripe     += fabs(diff1) * length;
            } // for emb

          }

          const uint64_t idx = (stripe-start_idx)*n_samples_r;
          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          //TFloat *dm_stripe = dm_stripes[stripe];

          // keep all writes in a single place, to maximize GPU warp performance
          dm_stripe[k] += my_stripe;
       } 

      } // for ik
     } // for stripe
    } // for sk

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

    const unsigned int step_size = su::UnifracVawUnnormalizedWeightedTask<TFloat>::step_size;
    const unsigned int sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // point of thread
#ifdef _OPENACC
    const unsigned int acc_vector_size = su::UnifracVawUnnormalizedWeightedTask<TFloat>::acc_vector_size;
#pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,lengths) async
#endif
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

    bool * const __restrict__ zcheck = this->zcheck;
    TFloat * const __restrict__ sums = this->sums;

    const unsigned int step_size = su::UnifracNormalizedWeightedTask<TFloat>::step_size;
    const unsigned int sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // check for zero values and pre-compute single column sums
#ifdef _OPENACC
#pragma acc parallel loop present(embedded_proportions,lengths,zcheck,sums)
#endif
    for(unsigned int k=0; k<n_samples; k++) {
            bool all_zeros=true;
            TFloat my_sum = 0.0;

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat u1 = embedded_proportions[offset + k];
                my_sum += u1*lengths[emb];
                all_zeros = all_zeros && (u1==0.0);
            }

            sums[k]     = my_sum;
            zcheck[k] = all_zeros;
    }

    // point of thread
#ifdef _OPENACC
    const unsigned int acc_vector_size = su::UnifracNormalizedWeightedTask<TFloat>::acc_vector_size;
#pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths,zcheck,sums) async
#endif
    for(unsigned int sk = 0; sk < sample_steps ; sk++) {
     for(unsigned int stripe = start_idx; stripe < stop_idx; stripe++) {
      for(unsigned int ik = 0; ik < step_size ; ik++) {
       const unsigned int k = sk*step_size + ik;

       if (k>=n_samples) continue; // past the limit

       const bool zcheck_k = zcheck[sk]; // due to loop collapse in ACC, must load in here

       const unsigned int l1 = (k + stripe + 1)%n_samples; // wraparound

       const bool allzero_k = zcheck[k];
       const bool allzero_l1 = zcheck[l1];

       if (allzero_k && allzero_l1) {
         // nothing to do, would have to add 0
       } else {
          const uint64_t idx = (stripe-start_idx) * n_samples_r;
          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
          //TFloat *dm_stripe = dm_stripes[stripe];
          //TFloat *dm_stripe_total = dm_stripes_total[stripe];

          // the totals can always use the distributed property
          dm_stripe_total[k] += sums[k] + sums[l1];
   
          TFloat my_stripe;

          if (allzero_k || allzero_l1) {
            // one side has all zeros
            // we can use the distributed property, and use the pre-computed values

            const unsigned int ridx = (allzero_k) ? l1 : // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
                                                    k;   // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0

            // keep reads in the same place to maximize GPU warp performance
            my_stripe = sums[ridx];
          
          } else {
            // both sides non zero, use the explicit but slow approach

            my_stripe = 0.0;

#pragma acc loop seq
            for (unsigned int emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat u1 = embedded_proportions[offset + k];
                TFloat v1 = embedded_proportions[offset + l1];
                TFloat diff1 = u1 - v1;
                TFloat length = lengths[emb];

                my_stripe     += fabs(diff1) * length;
            }

          }

          // keep all writes in a single place, to maximize GPU warp performance
          dm_stripe[k]       += my_stripe;
       }

      } // for ik
     } // for stripe
    } // for sk

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

    const unsigned int step_size = su::UnifracVawNormalizedWeightedTask<TFloat>::step_size;
    const unsigned int sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // point of thread
#ifdef _OPENACC
    const unsigned int acc_vector_size = su::UnifracVawNormalizedWeightedTask<TFloat>::acc_vector_size;
#pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
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

    const unsigned int step_size = su::UnifracGeneralizedTask<TFloat>::step_size;
    const unsigned int sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // point of thread
#ifdef _OPENACC
    const unsigned int acc_vector_size = su::UnifracGeneralizedTask<TFloat>::acc_vector_size;
#pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths) async
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

    const unsigned int step_size = su::UnifracVawGeneralizedTask<TFloat>::step_size;
    const unsigned int sample_steps = (n_samples+(step_size-1))/step_size; // round up
    // quick hack, to be finished

    // point of thread
#ifdef _OPENACC
    const unsigned int acc_vector_size = su::UnifracVawGeneralizedTask<TFloat>::acc_vector_size;
#pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
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
    const uint64_t * const __restrict__ embedded_proportions = this->embedded_proportions;
    TFloat * const __restrict__ dm_stripes_buf = this->dm_stripes.buf;
    TFloat * const __restrict__ dm_stripes_total_buf = this->dm_stripes_total.buf;

    TFloat * const __restrict__ sums = this->sums;

    const unsigned int step_size = su::UnifracUnweightedTask<TFloat>::step_size;
    const unsigned int sample_steps = (n_samples+(step_size-1))/step_size; // round up

    const unsigned int filled_embs_els = filled_embs/64;
    const unsigned int filled_embs_rem = filled_embs%64; 

    const unsigned int filled_embs_els_round = (filled_embs+63)/64;


    // pre-compute sums of length elements, since they are likely to be accessed many times
    // We will use a 8-bit map, to keep it small enough to keep in L1 cache
#ifdef _OPENACC
#pragma acc parallel loop collapse(2) gang present(lengths,sums) async
#endif
    for (unsigned int emb_el=0; emb_el<filled_embs_els; emb_el++) {
       for (unsigned int sub8=0; sub8<8; sub8++) {
          const unsigned int emb8 = emb_el*8+sub8;
          TFloat * __restrict__ psum = &(sums[emb8<<8]);
          const TFloat * __restrict__ pl   = &(lengths[emb8*8]);

#pragma acc loop vector
          // compute all the combinations for this block (8-bits total)
          // psum[0] = 0.0   // +0*pl[0]+0*pl[1]+0*pl[2]+...
          // psum[1] = pl[0] // +0*pl[1]+0*pl[2]+...
          // psum[2] = pl[1] // +0*pl[0]+0*pl[2]+
          // psum[2] = pl[0] + pl[1]
          // ...
          // psum[255] = pl[1] +.. + pl[7] // + 0*pl[0]
          // psum[255] = pl[0] +pl[1] +.. + pl[7]
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
#ifdef _OPENACC
#pragma acc parallel loop gang present(lengths,sums) async
#endif
       for (unsigned int sub8=0; sub8<8; sub8++) {
          // we are summing we have enough buffer in sums
          const unsigned int emb8 = emb_el*8+sub8;
          TFloat * __restrict__ psum = &(sums[emb8<<8]);

#pragma acc loop vector
          // compute all the combinations for this block, set to 0 any past the limit
          // as above
          for (unsigned int b8_i=0; b8_i<0x100; b8_i++) {
             TFloat val= 0;
             for (unsigned int li=(emb8*8); li<filled_embs; li++) {
               val += ((b8_i >>  (li-(emb8*8))) & 1) * lengths[li];
             }
             psum[b8_i] = val;
          }
        }
    }

    // point of thread
#ifdef _OPENACC
#pragma acc wait
    const unsigned int acc_vector_size = su::UnifracUnweightedTask<TFloat>::acc_vector_size;
#pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,sums) async
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

            bool did_update = false;
            TFloat my_stripe = 0.0;
            TFloat my_stripe_total = 0.0;

#pragma acc loop seq
            for (unsigned int emb_el=0; emb_el<filled_embs_els_round; emb_el++) {
                const uint64_t offset = n_samples_r * emb_el;
                const TFloat * __restrict__ psum = &(sums[emb_el*0x800]);

                uint64_t u1 = embedded_proportions[offset + k];
                uint64_t v1 = embedded_proportions[offset + l1];
                uint64_t o1 = u1 | v1;

                if (o1!=0) {  // zeros are prevalent
                    did_update=true;
                    uint64_t x1 = u1 ^ v1;

                    // Use the pre-computed sums
                    // Each range of 8 lengths has already been pre-computed and stored in psum
                    // Since embedded_proportions packed format is in 64-bit format for performance reasons
                    //    we need to add the 8 sums using the four 8-bits for addressing inside psum

                    my_stripe       += psum[              (x1 & 0xff)] + 
                                   psum[0x100+((x1 >>  8) & 0xff)] +
                                   psum[0x200+((x1 >> 16) & 0xff)] +
                                   psum[0x300+((x1 >> 24) & 0xff)] +
                                   psum[0x400+((x1 >> 32) & 0xff)] +
                                   psum[0x500+((x1 >> 40) & 0xff)] +
                                   psum[0x600+((x1 >> 48) & 0xff)] +
                                   psum[0x700+((x1 >> 56)       )];
                    my_stripe_total += psum[              (o1 & 0xff)] +
                                   psum[0x100+((o1 >>  8) & 0xff)] +
                                   psum[0x200+((o1 >> 16) & 0xff)] +
                                   psum[0x300+((o1 >> 24) & 0xff)] +
                                   psum[0x400+((o1 >> 32) & 0xff)] +
                                   psum[0x500+((o1 >> 40) & 0xff)] +
                                   psum[0x600+((o1 >> 48) & 0xff)] +
                                   psum[0x700+((o1 >> 56)       )];
                }
            }

            if (did_update) {
              dm_stripe[k]       += my_stripe;
              dm_stripe_total[k] += my_stripe_total;
            }

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

    const unsigned int step_size = su::UnifracVawUnweightedTask<TFloat>::step_size;
    const unsigned int sample_steps = (n_samples+(step_size-1))/step_size; // round up

    const unsigned int filled_embs_els = (filled_embs+31)/32; // round up

    // point of thread
#ifdef _OPENACC
    const unsigned int acc_vector_size = su::UnifracVawUnweightedTask<TFloat>::acc_vector_size;
#pragma acc parallel loop collapse(3) vector_length(acc_vector_size) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
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

