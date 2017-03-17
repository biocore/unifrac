#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include <unordered_map>
#include <stdlib.h>
#include <math.h>
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
        vec = (double*)malloc(sizeof(double) * defaultsize);
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

double** su::deconvolute_stripes_transpose(std::vector<double*> &stripes, uint32_t n) {
    // would be better to just do striped_to_condensed_form
    double **dm;
    dm = (double**)malloc(sizeof(double*) * n);
    for(unsigned int i = 0; i < n; i++) {
        dm[i] = (double*)malloc(sizeof(double) * n);
        dm[i][i] = 0;
    }

    for(unsigned int i = 0; i < n; i++) {
        double *vec = stripes[i];
        unsigned int k = 0;
        for(unsigned int col = i + 1; col < n; col++) {
            if(col < n) {
                dm[i][col] = vec[col];
                dm[col][i] = vec[col];
            } //else {
              //  dm[col][i] = vec[k];
              //  dm[row][col % n] = vec[k];
            //}
            //k++;
        }
    }
    return dm;
}

void _unnormalized_weighted_unifrac_transpose_task(std::vector<double*> &dm_stripes, 
                                         std::vector<double*> &dm_stripes_total,
                                         double* embedded_proportions,
                                         double length,
                                         uint32_t n_samples,
                                         unsigned int start,
                                         unsigned int stop) {
    double *dm_stripe;
    /* if striptes are transposed, then
     * for(unsigned int i = start; i < stop; i++) {
     *     u = proportions[i];
     *     stride_transpose = strides[i];
     *     for(unsigned int j = start + 1; j < n_samples; j++) {
     *         v = proportions[j];
     *         stride_transpose[j] += (u ^ v) * length
     */

    for(unsigned int stripe=start; stripe < stop; stripe++) {
        double u = embedded_proportions[stripe];
        double *st = dm_stripes[stripe];

        for(unsigned int j = stripe + 1; j < n_samples; j++) {
            double v = embedded_proportions[j];
            st[j] += fabs(u - v) * length;
        }
    }
}
void _unnormalized_weighted_unifrac_task(std::vector<double*> &dm_stripes, 
                                         std::vector<double*> &dm_stripes_total,
                                         double* embedded_proportions,
                                         double length,
                                         uint32_t n_samples,
                                         unsigned int start,
                                         unsigned int stop) {
    double *dm_stripe;
    for(unsigned int stripe=start; stripe < stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        for(unsigned int j = 0; j < n_samples; j++) {
            double u = embedded_proportions[j];
            double v = embedded_proportions[j + stripe + 1];
                
            dm_stripe[j] += fabs(u - v) * length;
        }
    }
}

void _normalized_weighted_unifrac_task(std::vector<double*> &dm_stripes, 
                                       std::vector<double*> &dm_stripes_total,
                                       double* embedded_proportions, 
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

        for(unsigned int j = 0; j < n_samples; j++) {
            double u = embedded_proportions[j];
            double v = embedded_proportions[j + stripe + 1];
               
            dm_stripe[j] += fabs(u - v) * length;
            dm_stripe_total[j] += fabs(u + v) * length;
        }
    }
}

void _unweighted_unifrac_task(std::vector<double*> &dm_stripes, 
                        std::vector<double*> &dm_stripes_total,
                        double* embedded_proportions, 
                        double length, 
                        uint32_t n_samples,
                        unsigned int start,
                        unsigned int stop) {
    double *dm_stripe;
    double *dm_stripe_total;
    
    // TODO: variable length stack allocated arrays are frowned upon so should
    // just malloc this
    bool bool_embedded[n_samples * 2];

    for(unsigned int i = 0; i < n_samples * 2; i++)
        bool_embedded[i] = embedded_proportions[i] > 0;

    for(unsigned int stripe = start; stripe < stop; stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < n_samples; j++) {
            bool u = bool_embedded[j];

            // can avoid embedding with
            // v = bool_embedded[(j + stripe + 1) % n_samples]
            // but modulus is expensive. can also be expressed as
            // v = bool_embedded[min(0, j + stripe + 1 - n_samples)]
            // ...or we just precompute an index such that v = bool_embedded[v_idx[j]]
            //
            bool v = bool_embedded[j + stripe + 1];
            
            dm_stripe[j] += (u ^ v) * length;
            dm_stripe_total[j] += (u | v) * length;
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
        case weighted_unnormalized_transpose:
            func = &_unnormalized_weighted_unifrac_transpose_task;
            break;
    }
    PropStack propstack(table.n_samples);

    uint32_t node;
    double *node_proportions;
    double **embedded_proportions = (double**)malloc(sizeof(double*) * 2);
    embedded_proportions[0] = (double*)malloc(sizeof(double) * table.n_samples * 2);
    embedded_proportions[1] = (double*)malloc(sizeof(double) * table.n_samples * 2);

    // the effect of this is to reserve (ceil(n_samples / 2))
    double length;
    //uint32_t n_rotations = end - start;

    // - 1 to avoid root   
    for(unsigned int k = 0; k < (tree.nparens / 2) - 1; k++) {
        node = tree.postorderselect(k);
        length = tree.lengths[node];

        // optional: embed this block into a prefetch thread. very easy. double buffer
        // already established for embedded_proportions
        node_proportions = propstack.pop(node);
        set_proportions(node_proportions, tree, node, table, propstack);
        embed_proportions(embedded_proportions[k % 2], node_proportions, table.n_samples);

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
        func(dm_stripes, dm_stripes_total, embedded_proportions[k % 2], length, 
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
    
    free(embedded_proportions[0]);
    free(embedded_proportions[1]);
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

    for(unsigned int i = 0; i < n_rotations; i++)
        dm_stripes[i] = (double*)calloc(sizeof(double), n_samples);
    
    return dm_stripes;
}

std::vector<double*> su::make_strides_transpose(unsigned int n_samples) {
    uint32_t n_rotations = (n_samples + 1) / 2;
    //std::vector<double*> dm_stripes(n_rotations);
    std::vector<double*> dm_stripes(n_samples);
    //for(unsigned int i = 0; i < n_rotations; i++)
    for(unsigned int i = 0; i < n_samples; i++)
        dm_stripes[i] = (double*)calloc(sizeof(double), n_samples);
    
    return dm_stripes;
}

