#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include <unordered_map>
#include <stdlib.h>

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

double* su::get_sample_counts(biom &table) {
    double *sample = (double*)malloc(sizeof(double) * table.n_obs);
    double *sample_counts = (double*)calloc(sizeof(double), table.n_samples);
    for(unsigned int i = 0; i < table.n_samples; i++) {
        table.get_sample_data(table.sample_ids[i], sample);
        for(unsigned int j = 0; j < table.n_obs; j++) {
            sample_counts[i] += sample[j];
        }
    }
    free(sample);
    return sample_counts;
}

void unnormalized_weighted_unifrac(std::vector<double*> &dm_stripes, 
                                   std::vector<double*> &dm_stripes_total,
                                   double* embedded_proportions, 
                                   double length, 
                                   uint32_t n_samples) {
    double *dm_stripe;

    //point of thread
    for(unsigned int stripe = 0; stripe < dm_stripes.size(); stripe++) {
        dm_stripe = dm_stripes[stripe];
        
        for(unsigned int j = 0; j < n_samples; j++) {
            double u = embedded_proportions[j];
            double v = embedded_proportions[j + stripe + 1];
                
            dm_stripe[j] += fabs(u - v) * length;
        }
    }
}

void normalized_weighted_unifrac(std::vector<double*> &dm_stripes, 
                                 std::vector<double*> &dm_stripes_total,
                                 double* embedded_proportions, 
                                 double length, 
                                 uint32_t n_samples) {
    double *dm_stripe;
    double *dm_stripe_total;

    // point of thread
    for(unsigned int stripe = 0; stripe < dm_stripes.size(); stripe++) {
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

void unweighted_unifrac(std::vector<double*> &dm_stripes, 
                        std::vector<double*> &dm_stripes_total,
                        double* embedded_proportions, 
                        double length, 
                        uint32_t n_samples) {
    double *dm_stripe;
    double *dm_stripe_total;
    bool bool_embedded[n_samples * 2];

    for(unsigned int i = 0; i < n_samples * 2; i++)
        bool_embedded[i] = embedded_proportions[i] > 0;

    // point of thread
    for(unsigned int stripe = 0; stripe < dm_stripes.size(); stripe++) {
        dm_stripe = dm_stripes[stripe];
        dm_stripe_total = dm_stripes_total[stripe];

        for(unsigned int j = 0; j < n_samples; j++) {
            bool u = bool_embedded[j];
            bool v = bool_embedded[j + stripe + 1];
               
            dm_stripe[j] += (u ^ v) * length;
            dm_stripe_total[j] += (u | v) * length;
        }
    }
}
// should return a DistanceMatrix...
double** su::unifrac(biom &table, BPTree &tree, Method unifrac_method) {
    void (*func)(std::vector<double*>&, std::vector<double*>&, double*, double, uint32_t);

    switch(unifrac_method) {
        case unweighted:
            func = &unweighted_unifrac;
            break;
        case weighted_normalized:
            func = &normalized_weighted_unifrac;
            break;
        case weighted_unnormalized:
            func = &unnormalized_weighted_unifrac;
            break;
    }
    
    PropStack propstack(table.n_samples);

    uint32_t node;
    double *node_proportions;
    double *embedded_proportions = (double*)malloc(sizeof(double) * table.n_samples * 2);
    double *sample_counts = get_sample_counts(table);

    // the effect of this is to reserve (ceil(n_samples / 2))
    uint32_t n_rotations = (table.n_samples + 1) / 2;
    std::vector<double*> dm_stripes(n_rotations);
    std::vector<double*> dm_stripes_total(n_rotations);
    double length;

    for(unsigned int i = 0; i < ((table.n_samples + 1) / 2); i++)
        dm_stripes[i] = (double*)calloc(sizeof(double), table.n_samples);
    
    if(unifrac_method == weighted_normalized || unifrac_method == unweighted) {
        for(unsigned int i = 0; i < ((table.n_samples + 1) / 2); i++) 
            dm_stripes_total[i] = (double*)calloc(sizeof(double), table.n_samples);
    } 

    for(unsigned int k = 0; k < (tree.nparens / 2); k++) {
        node = tree.postorderselect(k);
        length = tree.lengths[node];

        if(node == 0) // root
            break;

        // a prefetch can happen here in a separate thread while opencl or other
        // compute happens below
        node_proportions = propstack.pop(node);
        set_proportions(node_proportions, tree, node, table, propstack, sample_counts);
        embed_proportions(embedded_proportions, node_proportions, table.n_samples);
        // barrier for a prefetch

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

        // each stripe (or partial stripe) can be farmed out to opencl or other
        func(dm_stripes, dm_stripes_total, embedded_proportions, length, table.n_samples);
    }
    double **dm_unique;
    double **dm_total;
    dm_unique = deconvolute_stripes(dm_stripes, table.n_samples);
    
    if(unifrac_method == weighted_normalized || unifrac_method == unweighted) {
        dm_total = deconvolute_stripes(dm_stripes_total, table.n_samples);
        for(unsigned int row = 0; row < table.n_samples; row++) {
            for(unsigned int col = row + 1; col < table.n_samples; col++) {
                dm_unique[row][col] = dm_unique[row][col] / dm_total[row][col];
                dm_unique[col][row] = dm_unique[row][col];
            }
           free(dm_total[row]); 
        }
    }

    free(embedded_proportions);
    free(sample_counts);
    for(unsigned int i = 0; i < n_rotations; i++) {
        free(dm_stripes[i]);
        if(unifrac_method == weighted_normalized || unifrac_method == unweighted) 
            free(dm_stripes_total[i]);
    }

    return dm_unique; 
    /*
     * deconvolute the stripes
     */

    //return dm;
}

void su::set_proportions(double* props, 
                         BPTree &tree, 
                         uint32_t node, 
                         biom &table, 
                         PropStack &ps, 
                         double *sample_counts) {
    if(tree.isleaf(node)) {
       table.get_obs_data(tree.names[node], props);
       for(unsigned int i = 0; i < table.n_samples; i++)
           props[i] = props[i] / sample_counts[i];

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

