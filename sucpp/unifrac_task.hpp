#include "task_parameters.hpp"
#include <math.h>
#include <vector>
#include <stdint.h>
#include <stddef.h>

#ifndef __UNIFRAC_TASKS
#define __UNIFRAC_TASKS 1

namespace su {


#ifdef _OPENACC

  #ifndef SMALLGPU
  // defaultt on larger alignment, which improves performance on GPUs like V100
#define UNIFRAC_BLOCK 64
  #else
  // smaller GPUs prefer smaller allignment 
#define UNIFRAC_BLOCK 32
  #endif

#else

// CPUs don't need such a big alignment
#define UNIFRAC_BLOCK 8
#endif

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
    template<class TFloat>
    class UnifracTaskBase {
      public:
        UnifracTaskVector<TFloat> dm_stripes;
        UnifracTaskVector<TFloat> dm_stripes_total;

        const su::task_parameters* task_p;

        UnifracTaskBase(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const su::task_parameters* _task_p)
        : dm_stripes(_dm_stripes,_task_p), dm_stripes_total(_dm_stripes_total,_task_p), task_p(_task_p) {}

        // Note: not const, since they share a mutable state
        UnifracTaskBase(UnifracTaskBase &baseObj)
        : dm_stripes(baseObj.dm_stripes), dm_stripes_total(baseObj.dm_stripes_total), task_p(baseObj.task_p) {}

        virtual ~UnifracTaskBase() {}

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
     * task_p <task_parameters*> task specific parameters.
     */

    template<class TFloat>
    class UnifracTask : public UnifracTaskBase<TFloat> {
      protected:
#ifdef _OPENACC

        // The parallel nature of GPUs needs a largish step
  #ifndef SMALLGPU
        // default to larger step, which makes a big difference for bigger GPUs like V100
        static const unsigned int step_size = 32;
  #else
        // smaller GPUs prefer a slightly smaller step
        static const unsigned int step_size = 16;
  #endif
#else
        // The serial nature of CPU cores prefers a small step
        static const unsigned int step_size = 4;
#endif

      public:
        const TFloat * const embedded_proportions;
        const unsigned int max_embs;

        UnifracTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const TFloat * _embedded_proportions, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTaskBase<TFloat>(_dm_stripes, _dm_stripes_total, _task_p)
        , embedded_proportions(_embedded_proportions), max_embs(_max_embs) {}

        UnifracTask(UnifracTaskBase<TFloat> &baseObj, const TFloat * _embedded_proportions, unsigned int _max_embs)
        : UnifracTaskBase<TFloat>(baseObj)
        , embedded_proportions(_embedded_proportions), max_embs(_max_embs) {}

      

       virtual ~UnifracTask() {}

       virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) = 0;
    };


    template<class TFloat>
    class UnifracUnnormalizedWeightedTask : public UnifracTask<TFloat> {
      public:
        UnifracUnnormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const TFloat * _embedded_proportions, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat>(_dm_stripes,_dm_stripes_total,_embedded_proportions,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracNormalizedWeightedTask : public UnifracTask<TFloat> {
      public:
        UnifracNormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const TFloat * _embedded_proportions, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat>(_dm_stripes,_dm_stripes_total,_embedded_proportions,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracUnweightedTask : public UnifracTask<TFloat> {
      public:
        UnifracUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const TFloat * _embedded_proportions, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat>(_dm_stripes,_dm_stripes_total,_embedded_proportions,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracGeneralizedTask : public UnifracTask<TFloat> {
      public:
        UnifracGeneralizedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const TFloat * _embedded_proportions, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat>(_dm_stripes,_dm_stripes_total,_embedded_proportions,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };

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
     * task_p <task_parameters*> task specific parameters.
     */
    template<class TFloat>
    class UnifracVawTask : public UnifracTaskBase<TFloat> {
      protected:
#ifdef _OPENACC
        // The parallel nature of GPUs needs a largish step
  #ifndef SMALLGPU
        // default to larger step, which makes a big difference for bigger GPUs like V100
        static const unsigned int step_size = 32;
  #else
        // smaller GPUs prefer a slightly smaller step
        static const unsigned int step_size = 16;
  #endif
#else
        // The serial nature of CPU cores prefers a small step
        static const unsigned int step_size = 4;
#endif

      public:
        const TFloat * const embedded_proportions;
        const TFloat * const embedded_counts;
        const TFloat * const sample_total_counts;
        const unsigned int max_embs;

        UnifracVawTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const TFloat * _embedded_proportions, const TFloat * _embedded_counts, const TFloat * _sample_total_counts,
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTaskBase<TFloat>(_dm_stripes, _dm_stripes_total, _task_p)
        , embedded_proportions(_embedded_proportions), embedded_counts(_embedded_counts), sample_total_counts(_sample_total_counts), max_embs(_max_embs) {}

        UnifracVawTask(UnifracTaskBase<TFloat> &baseObj, 
                    const TFloat * _embedded_proportions, const TFloat * _embedded_counts, const TFloat * _sample_total_counts, unsigned int _max_embs)
        : UnifracTaskBase<TFloat>(baseObj)
        , embedded_proportions(_embedded_proportions), embedded_counts(_embedded_counts), sample_total_counts(_sample_total_counts), max_embs(_max_embs) {}



       virtual ~UnifracVawTask() {}

       virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) = 0;
    };

    template<class TFloat>
    class UnifracVawUnnormalizedWeightedTask : public UnifracVawTask<TFloat> {
      public:
        UnifracVawUnnormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const TFloat * _embedded_proportions, const TFloat * _embedded_counts, const TFloat * _sample_total_counts, 
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat>(_dm_stripes,_dm_stripes_total,_embedded_proportions,_embedded_counts,_sample_total_counts,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracVawNormalizedWeightedTask : public UnifracVawTask<TFloat> {
      public:
        UnifracVawNormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const TFloat * _embedded_proportions, const TFloat * _embedded_counts, const TFloat * _sample_total_counts, 
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat>(_dm_stripes,_dm_stripes_total,_embedded_proportions,_embedded_counts,_sample_total_counts,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracVawUnweightedTask : public UnifracVawTask<TFloat> {
      public:
        UnifracVawUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const TFloat * _embedded_proportions, const TFloat * _embedded_counts, const TFloat * _sample_total_counts, 
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat>(_dm_stripes,_dm_stripes_total,_embedded_proportions,_embedded_counts,_sample_total_counts,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };
    template<class TFloat>
    class UnifracVawGeneralizedTask : public UnifracVawTask<TFloat> {
      public:
        UnifracVawGeneralizedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total,
                    const TFloat * _embedded_proportions, const TFloat * _embedded_counts, const TFloat * _sample_total_counts, 
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat>(_dm_stripes,_dm_stripes_total,_embedded_proportions,_embedded_counts,_sample_total_counts,_max_embs,_task_p) {}

        virtual void run(unsigned int filled_embs, const TFloat * __restrict__ length) {_run(filled_embs, length);}

        void _run(unsigned int filled_embs, const TFloat * __restrict__ length);
    };

}

#endif
