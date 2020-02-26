#include <stack>
#include <vector>
#include <unordered_map>
#include <thread>
#include "unifrac_task.hpp"
#include <pthread.h>

#ifndef __UNIFRAC
    namespace su {
        enum Method {unweighted, weighted_normalized, weighted_unnormalized, generalized};
        
        class PropStack {
            private:
                std::stack<double*> prop_stack;
                std::unordered_map<uint32_t, double*> prop_map;
                uint32_t defaultsize;
            public:
                PropStack(uint32_t vecsize);
                ~PropStack();
                double* pop(uint32_t i);
                void push(uint32_t i);
                double* get(uint32_t i);
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
        void stripes_to_condensed_form(std::vector<double*> &stripes, uint32_t n, double* &cf, unsigned int start, unsigned int stop);
        void set_proportions(double* props, 
                             BPTree &tree, uint32_t node, 
                             biom &table, 
                             PropStack &ps,
                             bool normalize = true);
        std::vector<double*> make_strides(unsigned int n_samples);
        inline void embed_proportions(double* out, double* in, uint32_t n) {
#pragma acc parallel loop present(out) copyin(in[:n])
            for(unsigned int i = 0; i < n; i++) {
                double val = in[i];
                out[i] = val;
                out[i + n] = val;
            }
        }

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
