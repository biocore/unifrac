/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "task_parameters.hpp"
#include <math.h>
#include <vector>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>


#ifndef __UNIFRAC_TASKS
#define __UNIFRAC_TASKS 1

#ifdef _OPENACC

#define SUCMP_NM su_acc

  #ifndef SMALLGPU
  // defaultt on larger alignment, which improves performance on GPUs like V100
#define UNIFRAC_BLOCK 64
  #else
  // smaller GPUs prefer smaller allignment 
#define UNIFRAC_BLOCK 32
  #endif

#else

#define SUCMP_NM su_cpu


// CPUs don't need such a big alignment
#define UNIFRAC_BLOCK 16
#endif

namespace SUCMP_NM {

    // Note: This adds a copy, which is suboptimal
    //       But was the easiest way to get a contiguous buffer
    //       And it does allow for fp32 compute, when desired
    template<class TFloat>
    class UnifracTaskVector {
    private:
      std::vector<double*> &dm_stripes;
      const su::task_parameters* const task_p;

    public:
      const unsigned int start_idx;
      const unsigned int n_samples;
      const uint64_t  n_samples_r;
      TFloat* const buf;

      UnifracTaskVector(std::vector<double*> &_dm_stripes, const su::task_parameters* _task_p)
      : dm_stripes(_dm_stripes), task_p(_task_p)
      , start_idx(task_p->start), n_samples(task_p->n_samples)
      , n_samples_r(((n_samples + UNIFRAC_BLOCK-1)/UNIFRAC_BLOCK)*UNIFRAC_BLOCK) // round up
      , buf((dm_stripes[start_idx]==NULL) ? NULL : new TFloat[n_samples_r*(task_p->stop-start_idx)]) // dm_stripes could be null, in which case keep it null
      {
        TFloat* const ibuf = buf;
        if (ibuf != NULL) {
#ifdef _OPENACC
          const uint64_t bufels = n_samples_r * (task_p->stop-start_idx);
#endif
          for(unsigned int stripe=start_idx; stripe < task_p->stop; stripe++) {
             double * dm_stripe = dm_stripes[stripe];
             TFloat * buf_stripe = this->operator[](stripe);
             for(unsigned int j=0; j<n_samples; j++) {
                // Note: We could probably just initialize to zero
                buf_stripe[j] = dm_stripe[j];
             }
             for(unsigned int j=n_samples; j<n_samples_r; j++) {
                // Avoid NaNs
                buf_stripe[j] = 0.0;
             }
           }
#ifdef _OPENACC
#pragma acc enter data copyin(ibuf[:bufels])
#endif    
        }
      }

      TFloat * operator[](unsigned int idx) { return buf+((idx-start_idx)*n_samples_r);}
      const TFloat * operator[](unsigned int idx) const { return buf+((idx-start_idx)*n_samples_r);}


      ~UnifracTaskVector()
      {
        TFloat* const ibuf = buf;
        if (ibuf != NULL) {
#ifdef _OPENACC
          const uint64_t bufels = n_samples_r * (task_p->stop-start_idx); 
#pragma acc exit data copyout(ibuf[:bufels])
#endif    
          for(unsigned int stripe=start_idx; stripe < task_p->stop; stripe++) {
             double * dm_stripe = dm_stripes[stripe];
             TFloat * buf_stripe = this->operator[](stripe);
             for(unsigned int j=0; j<n_samples; j++) {
              dm_stripe[j] = buf_stripe[j];
             }
          }
          delete [] buf;
        }
      }

    private:
      UnifracTaskVector() = delete;
      UnifracTaskVector operator=(const UnifracTaskVector&other) const = delete;
    };

    // Base task class to be shared by all tasks
    template<class TFloat, class TEmb>
    class UnifracTaskBase {
      public:
        UnifracTaskVector<TFloat> dm_stripes;
        UnifracTaskVector<TFloat> dm_stripes_total;

        const su::task_parameters* task_p;

        const unsigned int max_embs;
        TEmb * embedded_proportions;
#ifdef _OPENACC
       protected:
        // alternate buffer only needed in async environments, like openacc
        TEmb * embedded_proportions_alt; // used as temp
       public:
#endif

        UnifracTaskBase(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : dm_stripes(_dm_stripes,_task_p), dm_stripes_total(_dm_stripes_total,_task_p), task_p(_task_p)
        , max_embs(_max_embs)
        , embedded_proportions(initialize_embedded(dm_stripes.n_samples_r,_max_embs))
#ifdef _OPENACC
        , embedded_proportions_alt(initialize_embedded(dm_stripes.n_samples_r,_max_embs)) 
#endif
        {}

        /* remove
        // Note: not const, since they share a mutable state
        UnifracTaskBase(UnifracTaskBase &baseObj)
        : dm_stripes(baseObj.dm_stripes), dm_stripes_total(baseObj.dm_stripes_total), task_p(baseObj.task_p) {}
        */

        virtual ~UnifracTaskBase()
        {
#ifdef _OPENACC
          const uint64_t  n_samples_r = dm_stripes.n_samples_r;
          const uint64_t bsize = n_samples_r * get_emb_els(max_embs);
#pragma acc exit data delete(embedded_proportions_alt[:bsize])
#pragma acc exit data delete(embedded_proportions[:bsize])
          free(embedded_proportions_alt);
#endif
          free(embedded_proportions);
        }

        void sync_embedded_proportions(unsigned int filled_embs)
        {
#ifdef _OPENACC
          const uint64_t  n_samples_r = dm_stripes.n_samples_r;
          const uint64_t bsize = n_samples_r * get_emb_els(filled_embs);
#pragma acc update device(embedded_proportions[:bsize])
#endif
        }

        static unsigned int get_emb_els(unsigned int max_embs);

        static TEmb *initialize_embedded(const uint64_t  n_samples_r, unsigned int max_embs) {
          uint64_t bsize = n_samples_r * get_emb_els(max_embs);

          TEmb* buf = NULL;
          int err = posix_memalign((void **)&buf, 4096, sizeof(TEmb) * bsize);
          if(buf == NULL || err != 0) {
            fprintf(stderr, "Failed to allocate %zd bytes, err %d; [%s]:%d\n",
                    sizeof(TEmb) * bsize, err, __FILE__, __LINE__);
             exit(EXIT_FAILURE);
          }
#pragma acc enter data create(buf[:bsize])
          return buf;
        }

        void embed_proportions_range(const TFloat* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb);
        void embed_proportions(const TFloat* __restrict__ in, unsigned int emb) {embed_proportions_range(in,0,dm_stripes.n_samples,emb);}

        

        //
        // ===== Internal, do not use directly =======
        //


        // Just copy from one buffer to another
        // May convert between fp formats in the process (if TOut!=double)
        template<class TOut> void embed_proportions_range_straight(TOut* __restrict__ out, const TFloat* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) const
        {
          const unsigned int n_samples  = dm_stripes.n_samples;
          const uint64_t n_samples_r  = dm_stripes.n_samples_r;
          const uint64_t offset = emb * n_samples_r;

          for(unsigned int i = start; i < end; i++) {
            out[offset + i] = in[i-start];
          }

          if (end==n_samples) {
            // avoid NaNs
            for(unsigned int i = n_samples; i < n_samples_r; i++) {
              out[offset + i] = 0.0;
            }
          }
        }

        // packed bool
        // Compute (in[:]>0) on each element, and store only the boolean bit.
        // The output values are stored in a multi-byte format, one bit per emb index,
        //    so it will likely take multiple passes to store all the values
        //
        // Note: assumes we are processing emb in increasing order, starting from 0
        template<class TOut> void embed_proportions_range_bool(TOut* __restrict__ out, const TFloat* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) const
        {
          const unsigned int n_packed = sizeof(TOut)*8;// e.g. 32 for unit32_t
          const unsigned int n_samples  = dm_stripes.n_samples;
          const uint64_t n_samples_r  = dm_stripes.n_samples_r;
          // The output values are stored in a multi-byte format, one bit per emb index
          // Compute the element to store the bit into, as well as whichbit in that element 
          unsigned int emb_block = emb/n_packed; // beginning of the element  block
          unsigned int emb_bit = emb%n_packed;   // bit inside the elements
          const uint64_t offset = emb_block * n_samples_r;

          if  (emb_bit==0) {
            // assign for emb_bit==0, so it clears the other bits
            // assumes we processing emb in increasing order, starting from 0
            for(unsigned int i = start; i < end; i++) {            
              out[offset + i] = (in[i-start] > 0);
            }

            if (end==n_samples) {
              // avoid NaNs
              for(unsigned int i = n_samples; i < n_samples_r; i++) {
                out[offset + i] = 0;
              }
            }
          } else {
            // just update my bit
            for(unsigned int i = start; i < end; i++) {
              out[offset + i] |= (TOut(in[i-start] > 0) << emb_bit);
            }

            // the rest of the els are already OK
          }
        }
    };

    // straight embeded_proportions
    template<> inline void UnifracTaskBase<double,double>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_straight(embedded_proportions,in,start,end,emb);}
    template<> inline void UnifracTaskBase<double,float>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_straight(embedded_proportions,in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,float>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_straight(embedded_proportions,in,start,end,emb);}
    template<> inline  unsigned int UnifracTaskBase<double,double>::get_emb_els(unsigned int max_embs) {return max_embs;}
    template<> inline  unsigned int UnifracTaskBase<double,float>::get_emb_els(unsigned int max_embs) {return max_embs;}
    template<> inline  unsigned int UnifracTaskBase<float,float>::get_emb_els(unsigned int max_embs) {return max_embs;}

    //packed bool embeded_proportions
    template<> inline void UnifracTaskBase<double,uint32_t>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_bool(embedded_proportions,in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,uint32_t>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_bool(embedded_proportions,in,start,end,emb);}
    template<> inline  unsigned int UnifracTaskBase<double,uint32_t>::get_emb_els(unsigned int max_embs) {return (max_embs+31)/32;}
    template<> inline  unsigned int UnifracTaskBase<float,uint32_t>::get_emb_els(unsigned int max_embs) {return (max_embs+31)/32;}

    template<> inline void UnifracTaskBase<double,uint64_t>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_bool(embedded_proportions,in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,uint64_t>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_bool(embedded_proportions,in,start,end,emb);}
    template<> inline  unsigned int UnifracTaskBase<double,uint64_t>::get_emb_els(unsigned int max_embs) {return (max_embs+63)/64;}
    template<> inline  unsigned int UnifracTaskBase<float,uint64_t>::get_emb_els(unsigned int max_embs) {return (max_embs+63)/64;}

    /* void unifrac tasks
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
     * task_p <task_parameters*> task specific parameters.
     */

    template<class TFloat, class TEmb>
    class UnifracTask : public UnifracTaskBase<TFloat,TEmb> {
      protected:
        // Use one cache line on CPU
        // On GPU, shaing a cache line is actually a good thing
        static const unsigned int step_size = 16*4/sizeof(TFloat);

#ifdef _OPENACC
        // Use as big vector size as we can, to maximize cache line reuse
        static const unsigned int acc_vector_size = 2048;
#endif

      public:

        UnifracTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTaskBase<TFloat,TEmb>(_dm_stripes, _dm_stripes_total, _max_embs, _task_p) {}

        /* delete
        UnifracTask(UnifracTaskBase<TFloat> &baseObj, const TEmb * _embedded_proportions, unsigned int _max_embs)
        : UnifracTaskBase<TFloat>(baseObj)
        , embedded_proportions(_embedded_proportions), max_embs(_max_embs) {}
        */
      

       virtual ~UnifracTask() {}

       virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) = 0;

      protected:
       static const unsigned int RECOMMENDED_MAX_EMBS_STRAIGHT = 128-16; // a little less to leave a bit of space of maxed-out L1
       // packed uses 32x less memory,so this should be 32x larger than straight... but there are additional structures, so use half of that
       static const unsigned int RECOMMENDED_MAX_EMBS_BOOL = 64*32;

    };


    template<class TFloat>
    class UnifracUnnormalizedWeightedTask : public UnifracTask<TFloat,TFloat> {
      public:
        static const unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<TFloat,TFloat>::RECOMMENDED_MAX_EMBS_STRAIGHT;

        UnifracUnnormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p)
        {
          const unsigned int n_samples = this->task_p->n_samples;

          zcheck = NULL;
          sums = NULL;
          posix_memalign((void **)&zcheck, 4096, sizeof(bool) * n_samples);
          posix_memalign((void **)&sums  , 4096, sizeof(TFloat) * n_samples);
#pragma acc enter data create(zcheck[:n_samples],sums[:n_samples])
        }

        virtual ~UnifracUnnormalizedWeightedTask()
        {
#ifdef _OPENACC
          const unsigned int n_samples = this->task_p->n_samples;
#pragma acc exit data delete(sums[:n_samples],zcheck[:n_samples])
#endif
          free(sums);
          free(zcheck);
        }

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
      protected:
        // temp buffers
        bool     *zcheck;
        TFloat   *sums;
    };
    template<class TFloat>
    class UnifracNormalizedWeightedTask : public UnifracTask<TFloat,TFloat> {
      public:
        static const unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<TFloat,TFloat>::RECOMMENDED_MAX_EMBS_STRAIGHT;

        UnifracNormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p)
        {
          const unsigned int n_samples = this->task_p->n_samples;

          zcheck = NULL;
          sums = NULL;
          posix_memalign((void **)&zcheck, 4096, sizeof(bool) * n_samples);
          posix_memalign((void **)&sums  , 4096, sizeof(TFloat) * n_samples);
#pragma acc enter data create(zcheck[:n_samples],sums[:n_samples])
        }

        virtual ~UnifracNormalizedWeightedTask()
        {
#ifdef _OPENACC
          const unsigned int n_samples = this->task_p->n_samples;
#pragma acc exit data delete(sums[:n_samples],zcheck[:n_samples])
#endif
          free(sums);
          free(zcheck);
        }

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
      protected:
        // temp buffers
        bool     *zcheck;
        TFloat   *sums;
    };
    template<class TFloat>
    class UnifracUnweightedTask : public UnifracTask<TFloat,uint64_t> {
      public:
        static const unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<TFloat,uint64_t>::RECOMMENDED_MAX_EMBS_BOOL;

        // Note: _max_emb MUST be multiple of 64
        UnifracUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat, uint64_t>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) 
        {
          const unsigned int bsize = _max_embs*(0x400/32);
          sums = NULL;
          posix_memalign((void **)&sums, 4096, sizeof(TFloat) * bsize);
#pragma acc enter data create(sums[:bsize])
        }

        virtual ~UnifracUnweightedTask()
        {
#ifdef _OPENACC
           const unsigned int bsize = this->max_embs*(0x400/32);
#pragma acc exit data delete(sums[:bsize])
#endif
          free(sums);
        }

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
      private:
        TFloat *sums; // temp buffer
    };
    template<class TFloat>
    class UnifracGeneralizedTask : public UnifracTask<TFloat,TFloat> {
      public:
        static const unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<TFloat,TFloat>::RECOMMENDED_MAX_EMBS_STRAIGHT;

        UnifracGeneralizedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };

    /* void unifrac_vaw tasks
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
     * task_p <task_parameters*> task specific parameters.
     */
    template<class TFloat, class TEmb>
    class UnifracVawTask : public UnifracTaskBase<TFloat,TEmb> {
      protected:
#ifdef _OPENACC
        // The parallel nature of GPUs needs a largish step
  #ifndef SMALLGPU
        // default to larger step, which makes a big difference for bigger GPUs like V100
        static const unsigned int step_size = 32;
        // keep the vector size just big enough to keep the used emb array inside the 32k buffer
        static const unsigned int acc_vector_size = 32*32*8/sizeof(TFloat);
  #else
        // smaller GPUs prefer a slightly smaller step
        static const unsigned int step_size = 16;
        // keep the vector size just big enough to keep the used emb array inside the 32k buffer
        static const unsigned int acc_vector_size = 32*32*8/sizeof(TFloat);
  #endif
#else
        // The serial nature of CPU cores prefers a small step
        static const unsigned int step_size = 4;
#endif

      public:
        TFloat * const embedded_counts;
        const TFloat * const sample_total_counts;

        static const unsigned int RECOMMENDED_MAX_EMBS = 128;

        UnifracVawTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const TFloat * _sample_total_counts,
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTaskBase<TFloat,TEmb>(_dm_stripes, _dm_stripes_total, _max_embs, _task_p)
        , embedded_counts(UnifracTaskBase<TFloat,TFloat>::initialize_embedded(this->dm_stripes.n_samples_r,_max_embs)), sample_total_counts(_sample_total_counts) {}


        /* delete
        UnifracVawTask(UnifracTaskBase<TFloat> &baseObj, 
                    const TEmb * _embedded_proportions, const TFloat * _sample_total_counts, unsigned int _max_embs)
        : UnifracTaskBase<TFloat>(baseObj)
        , embedded_proportions(_embedded_proportions), embedded_counts(initialize_embedded<TFloat>()), sample_total_counts(_sample_total_counts), max_embs(_max_embs) {}
        */


       virtual ~UnifracVawTask() {}

       void sync_embedded_counts(unsigned int filled_embs)
       {
#ifdef _OPENACC
          const uint64_t  n_samples_r = this->dm_stripes.n_samples_r;
          const uint64_t bsize = n_samples_r * filled_embs;
#pragma acc update device(embedded_counts[:bsize])
#endif
       }

       void sync_embedded(unsigned int filled_embs) { this->sync_embedded_proportions(filled_embs); this->sync_embedded_counts(filled_embs);}

        void embed_range(const TFloat* __restrict__ in_proportions, const TFloat* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb) {
          this->embed_proportions_range(in_proportions,start,end,emb);
          this->embed_proportions_range_straight(this->embedded_counts,in_counts,start,end,emb);
        }
        void embed(const TFloat* __restrict__ in_proportions, const double* __restrict__ in_counts, unsigned int emb) { embed_range(in_proportions,in_counts,0,this->dm_stripes.n_samples,emb);}

       virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) = 0;
    };

    template<class TFloat>
    class UnifracVawUnnormalizedWeightedTask : public UnifracVawTask<TFloat,TFloat> {
      public:
        UnifracVawUnnormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const TFloat * _sample_total_counts, 
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_sample_total_counts,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracVawNormalizedWeightedTask : public UnifracVawTask<TFloat,TFloat> {
      public:
        UnifracVawNormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const TFloat * _sample_total_counts, 
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_sample_total_counts,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracVawUnweightedTask : public UnifracVawTask<TFloat,uint32_t> {
      public:
        UnifracVawUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const TFloat * _sample_total_counts, 
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,uint32_t>(_dm_stripes,_dm_stripes_total,_sample_total_counts,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracVawGeneralizedTask : public UnifracVawTask<TFloat,TFloat> {
      public:
        UnifracVawGeneralizedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total,
                    const TFloat * _sample_total_counts, 
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_sample_total_counts,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };

}

#endif
