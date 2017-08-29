#include <stdint.h>

#ifndef __su_task_parameters
    #ifdef __cplusplus
    namespace su {
    #endif

        /* task specific compute parameters
         *
         * n_samples <int> the number of samples being processed
         * start <uint> the first stride to process
         * stop <uint> the last stride to process
         * tid <uint> the thread identifier
         * g_unifrac_alpha <double> an alpha value for generalized unifrac
         */
        struct task_parameters {
           uint32_t n_samples;          // number of samples
           unsigned int start;          // starting stripe
           unsigned int stop;           // stopping stripe
           unsigned int tid;            // thread ID
           
           // task specific arguments below
           double g_unifrac_alpha;      // generalized unifrac alpha
        };
    
    #ifdef __cplusplus
    }
    #endif

#define __su_task_parameters
#endif

