#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include <unordered_map>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <thread>

#if defined(__linux__)
    #include <sched.h>
    #include <pthread.h>
    // block below from http://bytefreaks.net/programming-2/cc-set-affinity-to-threads-example-code
    // note: comments adapted from the noted URL
    #include <errno.h>
     
    // The <errno.h> header file defines the integer variable errno, which is set by system calls and some library functions in the event of an error to indicate what went wrong.
    #define print_error_then_terminate(en, msg) \
      do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)
    #define print_perror_then_terminate(msg) \
      do { perror(msg); exit(EXIT_FAILURE); } while (0)
     
      struct thread_info {
     
        pthread_t thread_id; // ID returned by pthread_create()
        int core_id; // Core ID we want this pthread to set its affinity to
      };
 
#define SUCCESS_MSG "Successfully set thread %lu to affinity to CPU %d\n"
#define FAILURE_MSG "Failed to set thread %lu to affinity to CPU %d\n"
// end block
#endif


using namespace su;


PropStack::PropStack(uint32_t vecsize) {
    defaultsize = vecsize;
    prop_stack = std::stack<double*>();
    prop_map = std::unordered_map<uint32_t, double*>();

    prop_map.reserve(1000);
}

PropStack::~PropStack() {
    double *vec;
    // drain stack
    for(unsigned int i = 0; i < prop_stack.size(); i++) {
        vec = prop_stack.top();
        prop_stack.pop();
        free(vec);
    }
    
    // drain the map
    for(auto it = prop_map.begin(); it != prop_map.end(); it++) {
        vec = it->second;
        free(vec);
    }
    prop_map.clear();
}

double* PropStack::get(uint32_t i) {
    return prop_map[i];
}

void PropStack::push(uint32_t node) {
    double* vec = prop_map[node];
    prop_map.erase(node);
    prop_stack.push(vec);
}

double* PropStack::pop(uint32_t node) {
    /*
     * if we don't have any available vectors, create one
     * add it to our record of known vectors so we can track our mallocs
     */
    double *vec;
    if(prop_stack.empty()) {
        posix_memalign((void **)&vec, 32, sizeof(double) * defaultsize);
    }
    else {
        vec = prop_stack.top();
        prop_stack.pop();
    }

    prop_map[node] = vec;
    return vec;
}

double** su::deconvolute_stripes(std::vector<double*> &stripes, uint32_t n) {
    // would be better to just do striped_to_condensed_form
    double **dm;
    dm = (double**)malloc(sizeof(double*) * n);
    for(unsigned int i = 0; i < n; i++) {
        dm[i] = (double*)malloc(sizeof(double) * n);
        dm[i][i] = 0;
    }

    for(unsigned int i = 0; i < stripes.size(); i++) {
        double *vec = stripes[i];
        unsigned int k = 0;
        for(unsigned int row = 0, col = i + 1; row < n; row++, col++) {
            if(col < n) {
                dm[row][col] = vec[k];
                dm[col][row] = vec[k];
            } else {
                dm[col % n][row] = vec[k];
                dm[row][col % n] = vec[k];
            }
            k++;
        }
    }
    return dm;
}

void _unnormalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                         std::vector<double*> &__restrict__ dm_stripes_total,
                                         double* __restrict__ embedded_proportions,
                                         double length,
                                         uint32_t n_samples,
                                         unsigned int start,
                                         unsigned int stop) {
    double *dm_stripe;
    int j;
    //__m256d ymm0, ymm1; 
    for(int stripe=start; stripe < stop; stripe++) {
        dm_stripe = dm_stripes[stripe];

        /* intrinsics yield about a 2x reduction in runtime on llvm. they
         * were not effective on linux gcc 4.9.1 or 4.9.2. it is unclear
         * if they would be effective on other versions of gcc.
         *
         * one reason they help is that these for loops are not easily
         * autovectorizable. using the intrinsics effectively gets around
         * this. ...although, it also appears that loop unrolling works.
         *
         * it may make sense to revisit the inclusion of intriniscs, however
         * support must be tested at compile time, so it's rather annoying
         * at the moment. basically, we can't assume the presence of avx2.
         */
        //for(j = 0; j < n_samples / 4; j++) {
        //    int offset = j * 4;
        //    ymm0 = _mm256_sub_pd(_mm256_load_pd(embedded_proportions + offset),
        //                         _mm256_load_pd(embedded_proportions + (offset + stripe + 1))); // u - v
        //    ymm1 = _mm256_sub_pd(_mm256_set1_pd(0.0),
        //                         ymm0); // 0.0 - (u - v)
        //    ymm0 = _mm256_max_pd(ymm0, ymm1);  // abs(u - v)
        //    ymm0 = _mm256_fmadd_pd(ymm0,
        //                            _mm256_set1_pd(length),
        //                           _mm256_load_pd(dm_stripe + offset));  // fused mutiply add, dm_stripe[offset] += (abs(u - v) * length)
        //    _mm256_store_pd(dm_stripe + offset, ymm0);  //store
        //}

        //if((n_samples % 4) != 0) {
        //    for(int k = n_samples - (n_samples % 4); k < n_samples; k++) {
        //        double u = embedded_proportions[k];
        //        ////double v = alt_proportions[j];
        //        double v = embedded_proportions[k + stripe + 1];
        //        dm_stripe[k] += fabs(u - v) * length;
        //    }
        //}
        
        //for(int k = 0; k < n_samples; k++) {
        //    double u = embedded_proportions[k];
        //    double v = embedded_proportions[k + stripe + 1];
        //    dm_stripe[k] += fabs(u - v) * length;
        //}
        for(j = 0; j < n_samples / 4; j++) {
            int k = j * 4;
            double u1 = embedded_proportions[k];
            double u2 = embedded_proportions[k + 1];
            double u3 = embedded_proportions[k + 2];
            double u4 = embedded_proportions[k + 3];
            double v1 = embedded_proportions[k + stripe + 1];
            double v2 = embedded_proportions[k + stripe + 2];
            double v3 = embedded_proportions[k + stripe + 3];
            double v4 = embedded_proportions[k + stripe + 4];
            dm_stripe[k] += fabs(u1 - v1) * length;
            dm_stripe[k + 1] += fabs(u2 - v2) * length;
            dm_stripe[k + 2] += fabs(u3 - v3) * length;
            dm_stripe[k + 3] += fabs(u4 - v4) * length;
        }

        if((n_samples % 4) != 0) {
            for(int k = n_samples - (n_samples % 4); k < n_samples; k++) {
                double u = embedded_proportions[k];
                double v = embedded_proportions[k + stripe + 1];
                dm_stripe[k] += fabs(u - v) * length;
            }
        }
    }
}

void _normalized_weighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                                       std::vector<double*> &__restrict__ dm_stripes_total,
                                       double* __restrict__ embedded_proportions, 
                                       double length, 
                                       uint32_t n_samples,
                                       unsigned int start,
                                       unsigned int stop) {
    double *dm_stripe;
    double *dm_stripe_total;

    // point of thread
    for(unsigned int stripe = start; stripe < stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < n_samples / 4; j++) {
            int k = j * 4;
            double u1 = embedded_proportions[k];
            double u2 = embedded_proportions[k + 1];
            double u3 = embedded_proportions[k + 2];
            double u4 = embedded_proportions[k + 3];
         
            double v1 = embedded_proportions[k + stripe + 1];
            double v2 = embedded_proportions[k + stripe + 2];
            double v3 = embedded_proportions[k + stripe + 3];
            double v4 = embedded_proportions[k + stripe + 4];
               
            dm_stripe[k] += fabs(u1 - v1) * length;
            dm_stripe[k + 1] += fabs(u2 - v2) * length;
            dm_stripe[k + 2] += fabs(u3 - v3) * length;
            dm_stripe[k + 3] += fabs(u4 - v4) * length;

            dm_stripe_total[k] += fabs(u1 + v1) * length;
            dm_stripe_total[k + 1] += fabs(u2 + v2) * length;
            dm_stripe_total[k + 2] += fabs(u3 + v3) * length;
            dm_stripe_total[k + 3] += fabs(u4 + v4) * length;
        }

        if((n_samples % 4) != 0) {
            for(int k = n_samples - (n_samples % 4); k < n_samples; k++) {
                double u = embedded_proportions[k];
                double v = embedded_proportions[k + stripe + 1];
                   
                dm_stripe[k] += fabs(u - v) * length;
                dm_stripe_total[k] += fabs(u + v) * length;
            }
        }
    }
}

void _unweighted_unifrac_task(std::vector<double*> &__restrict__ dm_stripes, 
                              std::vector<double*> &__restrict__ dm_stripes_total,
                              double* __restrict__ embedded_proportions, 
                              double length,  
                              uint32_t n_samples,
                              unsigned int start,
                              unsigned int stop) {
    double *dm_stripe;
    double *dm_stripe_total;
    
    for(int stripe = start; stripe < stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(int j = 0; j < n_samples / 4; j++) {
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
        
        if((n_samples % 4) != 0) {
            for(int k = n_samples - (n_samples % 4); k < n_samples; k++) {
                int32_t u = embedded_proportions[k] > 0;
                int32_t v = embedded_proportions[k + stripe + 1] > 0;

                dm_stripe[k] += (u ^ v) * length;
                dm_stripe_total[k] += (u | v) * length;
            }
        }
    }
}


void progressbar(float progress) {
	// from http://stackoverflow.com/a/14539953
    //
    // could encapsulate into a classs for displaying time elapsed etc
	int barWidth = 70;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

void su::unifrac(biom &table,
                 BPTree &tree, 
                 Method unifrac_method,
                 std::vector<double*> &dm_stripes,
                 std::vector<double*> &dm_stripes_total,
                 unsigned int start, 
                 unsigned int end, 
                 unsigned int tid) {

    // bind self to a cpu (http://bytefreaks.net/programming-2/cc-set-affinity-to-threads-example-code)
    // note: comments adapted from the noted URL
#if defined(__linux__)
    const pthread_t pid = pthread_self();
    
    cpu_set_t cpuset; // cpu_set_t: This data set is a bitset where each bit represents a CPU.
    CPU_ZERO(&cpuset); // CPU_ZERO: This macro initializes the CPU set set to be the empty set.
    CPU_SET(tid, &cpuset); // CPU_SET: This macro adds cpu to the CPU set set.
    
    /* pthread_setaffinity_np: The pthread_setaffinity_np() function sets the CPU affinity 
     * mask of the thread thread to the CPU set pointed to by cpuset. If the call is successful, 
     * and the thread is not currently running on one of the CPUs in cpuset, 
     * then it is migrated to one of those CPUs.
     */
    const int set_result = pthread_setaffinity_np(pid, sizeof(cpu_set_t), &cpuset);
    if (set_result != 0) {
      print_error_then_terminate(set_result, "pthread_setaffinity_np");
    }
#endif

    void (*func)(std::vector<double*>&,  // dm_stripes
                 std::vector<double*>&,  // dm_stripes_total
                 double*,                // embedded_proportions
                 double,                 // length
                 uint32_t,               // number of samples
                 unsigned int,           // stripe start
                 unsigned int);          // stripe stop

    switch(unifrac_method) {
        case unweighted:
            func = &_unweighted_unifrac_task;
            break;
        case weighted_normalized:
            func = &_normalized_weighted_unifrac_task;
            break;
        case weighted_unnormalized:
            func = &_unnormalized_weighted_unifrac_task;
            break;
    }
    PropStack propstack(table.n_samples);

    uint32_t node;
    double *node_proportions;
    double *embedded_proportions; 
    double length;
	posix_memalign((void **)&embedded_proportions, 32, sizeof(double) * table.n_samples * 2);

    // thread local memory
    for(int i = start; i < end; i++){
        posix_memalign((void **)&dm_stripes[i], 32, sizeof(double) * table.n_samples);
        for(int j = 0; j < table.n_samples; j++)
            dm_stripes[i][j] = 0.;

        if(unifrac_method == unweighted || unifrac_method == weighted_normalized) {
            posix_memalign((void **)&dm_stripes_total[i], 32, sizeof(double) * table.n_samples);
            for(int j = 0; j < table.n_samples; j++)
                dm_stripes_total[i][j] = 0.;
        }
    }

    // - 1 to avoid root   
    for(unsigned int k = 0; k < (tree.nparens / 2) - 1; k++) {
        node = tree.postorderselect(k);
        length = tree.lengths[node];

        // optional: embed this block into a prefetch thread. very easy. double buffer
        // already established for embedded_proportions
        node_proportions = propstack.pop(node);
        set_proportions(node_proportions, tree, node, table, propstack);
        embed_proportions(embedded_proportions, node_proportions, table.n_samples);

        /*
         * The values in the example vectors correspond to index positions of an 
         * element in the resulting distance matrix. So, in the example below, 
         * the following can be interpreted:
         *
         * [0 1 2]
         * [1 2 3]
         *
         * As comparing the sample for row 0 against the sample for col 1, the
         * sample for row 1 against the sample for col 2, the sample for row 2
         * against the sample for col 3.
         *
         * In other words, we're computing stripes of a distance matrix. In the
         * following example, we're computing over 6 samples requiring 3 
         * stripes.
         *
         * A; stripe == 0
         * [0 1 2 3 4 5]
         * [1 2 3 4 5 0]
         *
         * B; stripe == 1
         * [0 1 2 3 4 5]
         * [2 3 4 5 0 1]
         *
         * C; stripe == 2
         * [0 1 2 3 4 5]
         * [3 4 5 0 1 2]
         *
         * The stripes end up computing the following positions in the distance
         * matrix.
         *
         * x A B C x x
         * x x A B C x
         * x x x A B C
         * C x x x A B
         * B C x x x A
         * A B C x x x
         *
         * However, we store those stripes as vectors, ie
         * [ A A A A A A ]
         *
         * We end up performing N / 2 redundant calculations on the last stripe 
         * (see C) but that is small over large N.  
         */
        func(dm_stripes, dm_stripes_total, embedded_proportions, length, 
             table.n_samples, start, end);
        
        // should make this compile-time support
        //if((tid == 0) && ((k % 1000) == 0))
 	    //    progressbar((float)k / (float)(tree.nparens / 2));       
    }
    
    if(unifrac_method == weighted_normalized || unifrac_method == unweighted) {
        for(unsigned int i = start; i < end; i++) {
            for(unsigned int j = 0; j < table.n_samples; j++) {
                dm_stripes[i][j] = dm_stripes[i][j] / dm_stripes_total[i][j];
            }
        }
    }
    
    free(embedded_proportions);
}

void su::set_proportions(double* props, 
                         BPTree &tree, 
                         uint32_t node, 
                         biom &table, 
                         PropStack &ps) {
    if(tree.isleaf(node)) {
       table.get_obs_data(tree.names[node], props);
       for(unsigned int i = 0; i < table.n_samples; i++)
           props[i] = props[i] / table.sample_counts[i];

    } else {
        unsigned int current = tree.leftchild(node);
        unsigned int right = tree.rightchild(node);
        double *vec;
        
        for(unsigned int i = 0; i < table.n_samples; i++)
            props[i] = 0;

        while(current <= right && current != 0) {
            vec = ps.get(current);  // pull from prop map
            ps.push(current);  // remove from prop map, place back on stack

            for(unsigned int i = 0; i < table.n_samples; i++)
                props[i] = props[i] + vec[i];

            current = tree.rightsibling(current);
        }
    }    
}

std::vector<double*> su::make_strides(unsigned int n_samples) {
    uint32_t n_rotations = (n_samples + 1) / 2;
    std::vector<double*> dm_stripes(n_rotations);

    for(unsigned int i = 0; i < n_rotations; i++) {
        double* tmp;
        posix_memalign((void **)&tmp, 32, sizeof(double) * n_samples);
        for(int j = 0; j < n_samples; j++)
            tmp[j] = 0.0;
        dm_stripes[i] = tmp;
    }    
    return dm_stripes;
}
