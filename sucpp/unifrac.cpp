#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include <unordered_map>
#include <cstdlib>
#include <algorithm>
#include <thread>


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
        if(vec == NULL) {
            fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                    sizeof(double) * defaultsize, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
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
    if(dm == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(double*) * n, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    for(unsigned int i = 0; i < n; i++) {
        dm[i] = (double*)malloc(sizeof(double) * n);
        if(dm[i] == NULL) {
            fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                    sizeof(double) * n, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
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
    //__m256d ymm0, ymm1; 
    for(unsigned int stripe=start; stripe < stop; stripe++) {
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
        }

        if((n_samples % 4) != 0) {
            for(unsigned int k = n_samples - (n_samples % 4); k < n_samples; k++) {
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
            for(unsigned int k = n_samples - (n_samples % 4); k < n_samples; k++) {
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
    
    for(unsigned int stripe = start; stripe < stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < n_samples / 4; j++) {
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
            for(unsigned int k = n_samples - (n_samples % 4); k < n_samples; k++) {
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

void initialize_embedded(double*& prop, const su::task_parameters* task_p) {
	posix_memalign((void **)&prop, 32, sizeof(double) * task_p->n_samples * 2);
    if(prop == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(double) * task_p->n_samples * 2, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
}

void initialize_sample_counts(double*& counts, const su::task_parameters* task_p, biom &table) {
	posix_memalign((void **)&counts, 32, sizeof(double) * task_p->n_samples * 2);
    if(counts == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(double) * task_p->n_samples, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    for(unsigned int i = 0; i < table.n_samples; i++) {
        counts[i] = table.sample_counts[i];
        counts[i + table.n_samples] = table.sample_counts[i];
    }
}

void initialize_stripes(std::vector<double*> &dm_stripes, 
                        std::vector<double*> &dm_stripes_total, 
                        Method unifrac_method, 
                        const su::task_parameters* task_p) {
    for(unsigned int i = task_p->start; i < task_p->stop; i++){
        posix_memalign((void **)&dm_stripes[i], 32, sizeof(double) * task_p->n_samples);
        if(dm_stripes[i] == NULL) {
            fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                    sizeof(double) * task_p->n_samples, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        for(unsigned int j = 0; j < task_p->n_samples; j++)
            dm_stripes[i][j] = 0.;

        if(unifrac_method == unweighted || unifrac_method == weighted_normalized) {
            posix_memalign((void **)&dm_stripes_total[i], 32, sizeof(double) * task_p->n_samples);
            if(dm_stripes_total[i] == NULL) {
                fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                        sizeof(double) * task_p->n_samples, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            for(unsigned int j = 0; j < task_p->n_samples; j++)
                dm_stripes_total[i][j] = 0.;
        }
    }
}


void su::unifrac(biom &table,
                 BPTree &tree, 
                 Method unifrac_method,
                 std::vector<double*> &dm_stripes,
                 std::vector<double*> &dm_stripes_total,
                 const su::task_parameters* task_p) {

    if(table.n_samples != task_p->n_samples) {
        fprintf(stderr, "Task and table n_samples not equal\n");
        exit(EXIT_FAILURE);
    }

    void (*func)(std::vector<double*>&,  // dm_stripes
                 std::vector<double*>&,  // dm_stripes_total
                 double*,                // embedded_proportions
                 double,                 // length
                 const su::task_parameters*);

    switch(unifrac_method) {
        case unweighted:
            func = &su::_unweighted_unifrac_task;
            break;
        case weighted_normalized:
            func = &su::_normalized_weighted_unifrac_task;
            break;
        case weighted_unnormalized:
            func = &su::_unnormalized_weighted_unifrac_task;
            break;
        case generalized:
            func = &su::_generalized_unifrac_task;
            break;
    }
    PropStack propstack(table.n_samples);

    uint32_t node;
    double *node_proportions;
    double *embedded_proportions; 
    double length;

    initialize_embedded(embedded_proportions, task_p);
    initialize_stripes(std::ref(dm_stripes), std::ref(dm_stripes_total), unifrac_method, task_p);

    // - 1 to avoid root   
    for(unsigned int k = 0; k < (tree.nparens / 2) - 1; k++) {
        node = tree.postorderselect(k);
        length = tree.lengths[node];

        node_proportions = propstack.pop(node);
        set_proportions(node_proportions, tree, node, table, propstack);
        embed_proportions(embedded_proportions, node_proportions, task_p->n_samples);

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
        func(dm_stripes, dm_stripes_total, embedded_proportions, length, task_p);
        
        // should make this compile-time support
        //if((tid == 0) && ((k % 1000) == 0))
 	    //    progressbar((float)k / (float)(tree.nparens / 2));       
    }
    
    if(unifrac_method == weighted_normalized || unifrac_method == unweighted || unifrac_method == generalized) {
        for(unsigned int i = task_p->start; i < task_p->stop; i++) {
            for(unsigned int j = 0; j < task_p->n_samples; j++) {
                dm_stripes[i][j] = dm_stripes[i][j] / dm_stripes_total[i][j];
            }
        }
    }
    
    free(embedded_proportions);
}

void su::unifrac_vaw(biom &table,
                     BPTree &tree, 
                     Method unifrac_method,
                     std::vector<double*> &dm_stripes,
                     std::vector<double*> &dm_stripes_total,
                     const su::task_parameters* task_p) {

    if(table.n_samples != task_p->n_samples) {
        fprintf(stderr, "Task and table n_samples not equal\n");
        exit(EXIT_FAILURE);
    }

    void (*func)(std::vector<double*>&,  // dm_stripes
                 std::vector<double*>&,  // dm_stripes_total
                 double*,                // embedded_proportions
                 double*,                // embedded_counts
                 double*,                // sample total counts
                 double,                 // length
                 const su::task_parameters*);

    switch(unifrac_method) {
        case unweighted:
            func = &su::_vaw_unweighted_unifrac_task;
            break;
        case weighted_normalized:
            func = &su::_vaw_normalized_weighted_unifrac_task;
            break;
        case weighted_unnormalized:
            func = &su::_vaw_unnormalized_weighted_unifrac_task;
            break;
        case generalized:
            func = &su::_vaw_generalized_unifrac_task;
            break;
    }
    PropStack propstack(table.n_samples);
    PropStack countstack(table.n_samples);

    uint32_t node;
    double *node_proportions;
    double *node_counts;
    double *embedded_proportions; 
    double *embedded_counts; 
    double *sample_total_counts;
    double length;

    initialize_embedded(embedded_proportions, task_p);
    initialize_embedded(embedded_counts, task_p);
    initialize_sample_counts(sample_total_counts, task_p, table);
    initialize_stripes(std::ref(dm_stripes), std::ref(dm_stripes_total), unifrac_method, task_p);

    // - 1 to avoid root   
    for(unsigned int k = 0; k < (tree.nparens / 2) - 1; k++) {
        node = tree.postorderselect(k);
        length = tree.lengths[node];

        node_proportions = propstack.pop(node);
        node_counts = countstack.pop(node);

        set_proportions(node_proportions, tree, node, table, propstack);
        set_proportions(node_counts, tree, node, table, countstack, false);

        embed_proportions(embedded_proportions, node_proportions, task_p->n_samples);
        embed_proportions(embedded_counts, node_counts, task_p->n_samples);

        func(dm_stripes, dm_stripes_total, embedded_proportions, embedded_counts, sample_total_counts, length, task_p);
    }
    
    if(unifrac_method == weighted_normalized || unifrac_method == unweighted || unifrac_method == generalized) {
        for(unsigned int i = task_p->start; i < task_p->stop; i++) {
            for(unsigned int j = 0; j < task_p->n_samples; j++) {
                dm_stripes[i][j] = dm_stripes[i][j] / dm_stripes_total[i][j];
            }
        }
    }
    
    free(embedded_proportions);
    free(embedded_counts);
    free(sample_total_counts);
}
void su::set_proportions(double* props, 
                         BPTree &tree, 
                         uint32_t node, 
                         biom &table, 
                         PropStack &ps,
                         bool normalize) {
    if(tree.isleaf(node)) {
       table.get_obs_data(tree.names[node], props);
       for(unsigned int i = 0; i < table.n_samples; i++) {
           props[i] = props[i];
           if(normalize)
               props[i] /= table.sample_counts[i];
       }

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
        if(tmp == NULL) {
            fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                    sizeof(double) * n_samples, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        for(unsigned int j = 0; j < n_samples; j++)
            tmp[j] = 0.0;
        dm_stripes[i] = tmp;
    }    
    return dm_stripes;
}
