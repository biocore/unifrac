#include <stack>
#include <vector>
#include <unordered_map>
#include <thread>
#include "unifrac_task.hpp"
#include <pthread.h>

#ifndef __UNIFRAC
    namespace su {
        enum Method {unweighted, weighted_normalized, weighted_unnormalized, generalized, unweighted_fp32, weighted_normalized_fp32, weighted_unnormalized_fp32, generalized_fp32};

        template<class TFloat> 
        class PropStack {
            private:
                std::stack<TFloat*> prop_stack;
                std::unordered_map<uint32_t, TFloat*> prop_map;
                uint32_t defaultsize;
            public:
                PropStack(uint32_t vecsize);
                ~PropStack();
                TFloat* pop(uint32_t i);
                void push(uint32_t i);
                TFloat* get(uint32_t i);
        };

        void faith_pd(biom &table, BPTree &tree, double* result);

        std::string test_table_ids_are_subset_of_tree(biom &table, BPTree &tree);
        void unifrac(biom &table, 
                     BPTree &tree, 
                     Method unifrac_method,
                     std::vector<double*> &dm_stripes,
                     std::vector<double*> &dm_stripes_total,
                     const task_parameters* task_p);
        
        void unifrac_vaw(biom &table, 
                         BPTree &tree, 
                         Method unifrac_method,
                         std::vector<double*> &dm_stripes,
                         std::vector<double*> &dm_stripes_total,
                         const task_parameters* task_p);
        
        double** deconvolute_stripes(std::vector<double*> &stripes, uint32_t n);

        class ManagedStripes {
        public:
           virtual ~ManagedStripes() {}
           virtual const double *get_stripe(uint32_t stripe) const = 0;
           virtual void release_stripe(uint32_t stripe) const = 0;
        };

        class MemoryStripes : public ManagedStripes {
        private:
           const double  * const * stripes;  // just a pointer, not owned
        public:
           MemoryStripes(const double  * const * _stripes) : stripes(_stripes) {}
           MemoryStripes(std::vector<double*> &_stripes) : stripes(_stripes.data()) {}
           MemoryStripes(const std::vector<double*> &_stripes) : stripes(_stripes.data()) {}
           MemoryStripes(const std::vector<const double*> &_stripes) : stripes(_stripes.data()) {}
           MemoryStripes(std::vector<const double*> &_stripes) : stripes(_stripes.data()) {}

           virtual const double *get_stripe(uint32_t stripe) const {return stripes[stripe];}
           virtual void release_stripe(uint32_t stripe) const {};
        };


        void stripes_to_condensed_form(std::vector<double*> &stripes, uint32_t n, double* cf, unsigned int start, unsigned int stop);

        // tile_size==0 means memory optimized
        template<class TReal> void stripes_to_matrix_T(const ManagedStripes &stripes, const uint32_t n_samples, const uint32_t n_stripes, TReal*  __restrict__ buf2d, uint32_t tile_size=0);
        void stripes_to_matrix(const ManagedStripes &stripes, const uint32_t n_samples, const uint32_t n_stripes, double*  __restrict__ buf2d, uint32_t tile_size=0);
        void stripes_to_matrix_fp32(const ManagedStripes &stripes, const uint32_t n_samples, const uint32_t n_stripes, float*  __restrict__ buf2d, uint32_t tile_size=0);


        template<class TReal> void condensed_form_to_matrix_T(const double*  __restrict__ cf, const uint32_t n, TReal*  __restrict__ buf2d);
        void condensed_form_to_matrix(const double*  __restrict__ cf, const uint32_t n, double*  __restrict__ buf2d);
        void condensed_form_to_matrix_fp32(const double*  __restrict__ cf, const uint32_t n, float*  __restrict__ buf2d);

        template<class TFloat>
        void set_proportions(TFloat* __restrict__ props, 
                             const BPTree &tree, uint32_t node, 
                             const biom &table, 
                             PropStack<TFloat> &ps,
                             bool normalize = true);

        template<class TFloat>
        void set_proportions_range(TFloat* __restrict__ props,
                                   const BPTree &tree, uint32_t node,
                                   const biom &table,unsigned int start, unsigned int end,
                                   PropStack<TFloat> &ps,
                                   bool normalize = true);

        std::vector<double*> make_strides(unsigned int n_samples);

        inline uint64_t comb_2(uint64_t N) {
            // based off of _comb_int_long
            // https://github.com/scipy/scipy/blob/v0.19.1/scipy/special/_comb.pyx
            
            // Compute binom(N, k) for integers.
            //
            // we're disregarding overflow as that practically should not
            // happen unless the number of samples processed is in excess
            // of 4 billion 
            uint64_t val, j, M, nterms;
            uint64_t k = 2;

            M = N + 1;
            nterms = k < (N - k) ? k : N - k;

            val = 1;

            for(j = 1; j < nterms + 1; j++) {
                val *= M - j;
                val /= j;
            }
            return val;
        }

        // process the stripes described by tasks
        void process_stripes(biom &table, 
                             BPTree &tree_sheared, 
                             Method method,
                             bool variance_adjust,
                             std::vector<double*> &dm_stripes, 
                             std::vector<double*> &dm_stripes_total,
                             std::vector<std::thread> &threads,
                             std::vector<su::task_parameters> &tasks);
    }
#define __UNIFRAC 1
#endif
