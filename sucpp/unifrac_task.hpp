#include "task_parameters.hpp"
#include <math.h>
#include <vector>
#include <stdint.h>
#include <stddef.h>

namespace su {

    // Note: This adds a copy, which is suboptimal
    //       But was the easiest way to get a contiguous buffer
    //       Future improvement welcome
    class UnifracTaskVector {
    private:
      std::vector<double*> &dm_stripes;
      const su::task_parameters* const task_p;

    public:
      const unsigned int start_idx;
      const unsigned int n_samples;
      double* const buf;

      UnifracTaskVector(std::vector<double*> &_dm_stripes, const su::task_parameters* _task_p)
      : dm_stripes(_dm_stripes), task_p(_task_p)
      , start_idx(task_p->start), n_samples(task_p->n_samples)
      , buf((dm_stripes[start_idx]==NULL) ? NULL : new double[n_samples*(task_p->stop-start_idx)]) // dm_stripes could be null, in which case keep it null
      {
        double* const ibuf = buf;
        if (ibuf != NULL) {
#ifdef _OPENACC
          unsigned int bufels = n_samples*(task_p->stop-start_idx);
#endif
          for(unsigned int stripe=start_idx; stripe < task_p->stop; stripe++) {
             double * dm_stripe = dm_stripes[stripe];
             double * buf_stripe = this->operator[](stripe);
             for(unsigned int j=0; j<n_samples; j++) {
                // Note: We could probably just initialize to zero
                buf_stripe[j] = dm_stripe[j];
             }
           }
#ifdef _OPENACC
#pragma acc enter data copyin(ibuf[:bufels])
#endif    
        }
      }

      double * operator[](unsigned int idx) { return buf+((idx-start_idx)*n_samples);}
      const double * operator[](unsigned int idx) const { return buf+((idx-start_idx)*n_samples);}


      ~UnifracTaskVector()
      {
        double* const ibuf = buf;
        if (ibuf != NULL) {
#ifdef _OPENACC
          unsigned int bufels = n_samples*(task_p->stop-start_idx);
#pragma acc exit data copyout(ibuf[:bufels])
#endif    
          for(unsigned int stripe=start_idx; stripe < task_p->stop; stripe++) {
             double * dm_stripe = dm_stripes[stripe];
             double * buf_stripe = this->operator[](stripe);
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
    class UnifracTaskBase {
      public:
        UnifracTaskVector dm_stripes;
        UnifracTaskVector dm_stripes_total;

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

    class UnifracTask : public UnifracTaskBase {
      public:
        const double * const embedded_proportions;

        UnifracTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const double * _embedded_proportions, const su::task_parameters* _task_p)
        : UnifracTaskBase(_dm_stripes, _dm_stripes_total, _task_p)
        , embedded_proportions(_embedded_proportions) {}

        UnifracTask(UnifracTaskBase &baseObj, const double * _embedded_proportions)
        : UnifracTaskBase(baseObj)
        , embedded_proportions(_embedded_proportions) {}

      

       virtual ~UnifracTask() {}

       virtual void run(double length) = 0;
    };


    class UnifracUnnormalizedWeightedTask : public UnifracTask {
      public:
        UnifracUnnormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const double * _embedded_proportions, const su::task_parameters* _task_p)
        : UnifracTask(_dm_stripes,_dm_stripes_total,_embedded_proportions,_task_p) {}

        virtual void run(double length) {_run(length);}

        void _run(double length);
    };
    class UnifracNormalizedWeightedTask : public UnifracTask {
      public:
        UnifracNormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const double * _embedded_proportions, const su::task_parameters* _task_p)
        : UnifracTask(_dm_stripes,_dm_stripes_total,_embedded_proportions,_task_p) {}

        virtual void run(double length) {_run(length);}

        void _run(double length);
    };
    class UnifracUnweightedTask : public UnifracTask {
      public:
        UnifracUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const double * _embedded_proportions, const su::task_parameters* _task_p)
        : UnifracTask(_dm_stripes,_dm_stripes_total,_embedded_proportions,_task_p) {}

        virtual void run(double length) {_run(length);}

        void _run(double length);
    };
    class UnifracGeneralizedTask : public UnifracTask {
      public:
        UnifracGeneralizedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, const double * _embedded_proportions, const su::task_parameters* _task_p)
        : UnifracTask(_dm_stripes,_dm_stripes_total,_embedded_proportions,_task_p) {}

        virtual void run(double length) {_run(length);}

        void _run(double length);
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
    class UnifracVawTask : public UnifracTaskBase {
      public:
        const double * const embedded_proportions;
        const double * const embedded_counts;
        const double * const sample_total_counts;

        UnifracVawTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _embedded_proportions, const double * _embedded_counts, const double * _sample_total_counts,
                    const su::task_parameters* _task_p)
        : UnifracTaskBase(_dm_stripes, _dm_stripes_total, _task_p)
        , embedded_proportions(_embedded_proportions), embedded_counts(_embedded_counts), sample_total_counts(_sample_total_counts) {}

        UnifracVawTask(UnifracTaskBase &baseObj, 
                    const double * _embedded_proportions, const double * _embedded_counts, const double * _sample_total_counts)
        : UnifracTaskBase(baseObj)
        , embedded_proportions(_embedded_proportions), embedded_counts(_embedded_counts), sample_total_counts(_sample_total_counts) {}



       virtual ~UnifracVawTask() {}

       virtual void run(double length) = 0;
    };

    class UnifracVawUnnormalizedWeightedTask : public UnifracVawTask {
      public:
        UnifracVawUnnormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _embedded_proportions, const double * _embedded_counts, const double * _sample_total_counts, 
                    const su::task_parameters* _task_p)
        : UnifracVawTask(_dm_stripes,_dm_stripes_total,_embedded_proportions,_embedded_counts,_sample_total_counts,_task_p) {}

        virtual void run(double length) {_run(length);}

        void _run(double length);
    };
    class UnifracVawNormalizedWeightedTask : public UnifracVawTask {
      public:
        UnifracVawNormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _embedded_proportions, const double * _embedded_counts, const double * _sample_total_counts, 
                    const su::task_parameters* _task_p)
        : UnifracVawTask(_dm_stripes,_dm_stripes_total,_embedded_proportions,_embedded_counts,_sample_total_counts,_task_p) {}

        virtual void run(double length) {_run(length);}

        void _run(double length);
    };
    class UnifracVawUnweightedTask : public UnifracVawTask {
      public:
        UnifracVawUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _embedded_proportions, const double * _embedded_counts, const double * _sample_total_counts, 
                    const su::task_parameters* _task_p)
        : UnifracVawTask(_dm_stripes,_dm_stripes_total,_embedded_proportions,_embedded_counts,_sample_total_counts,_task_p) {}

        virtual void run(double length) {_run(length);}

        void _run(double length);
    };
    class UnifracVawGeneralizedTask : public UnifracVawTask {
      public:
        UnifracVawGeneralizedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total,
                    const double * _embedded_proportions, const double * _embedded_counts, const double * _sample_total_counts, 
                    const su::task_parameters* _task_p)
        : UnifracVawTask(_dm_stripes,_dm_stripes_total,_embedded_proportions,_embedded_counts,_sample_total_counts,_task_p) {}

        virtual void run(double length) {_run(length);}

        void _run(double length);
    };

}
