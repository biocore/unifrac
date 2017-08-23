#include "api.hpp"
#include <cstring>

using namespace su;
using namespace std;

// https://stackoverflow.com/a/19841704/19741
bool is_file_exists(const char *fileName) {
    std::ifstream infile(fileName);
        return infile.good();
}


void destroy_stripes(vector<double*> &dm_stripes, vector<double*> &dm_stripes_total, unsigned int n_samples) {
    unsigned int n_rotations = (n_samples + 1) / 2;
    for(unsigned int i = 0; i < n_rotations; i++) {
        free(dm_stripes[i]);
        if(dm_stripes_total[i] != NULL) {
            free(dm_stripes_total[i]);
        }
    }
}

void su::initialize_mat(mat* &result, biom &table) {
    result = (mat*)malloc(sizeof(mat));
    result->n_samples = table.n_samples;
    
    // // // // // // // condensed form represents full distance matrix, refactor is not complete
    result->cf_size = table.n_samples * table.n_samples;
    result->is_square = true;
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);

    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table.sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table.sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
    }
}

void su::destroy_mat(mat* &result) {
    // // // // condensed form is currently full distance matrix
    for(unsigned int i = 0; i < result->n_samples; i++) {
        free(result->sample_ids[i]);
        free(result->condensed_form[i]);
    };
    free(result->sample_ids);
    free(result->condensed_form);
    free(result);
}

su::compute_status su::one_off(const char* biom_filename, const char* tree_filename, 
                               const char* unifrac_method, bool variance_adjust, double alpha,
                               unsigned int nthreads, mat* &result) {
    if(!is_file_exists(biom_filename)) {
        fprintf(stderr, "%s does not exist.\n", biom_filename);
        return su::table_missing;
    }

    if(!is_file_exists(tree_filename)) {
        fprintf(stderr, "%s does not exist.\n", biom_filename);
        return tree_missing;
    }

    Method method;
    if(std::strcmp(unifrac_method, "unweighted") == 0)
        method = su::unweighted;
    else if(std::strcmp(unifrac_method, "weighted_normalized") == 0)
        method = su::weighted_normalized;
    else if(std::strcmp(unifrac_method, "weighted_unnormalized") == 0)
        method = su::weighted_unnormalized;
    else if(std::strcmp(unifrac_method, "generalized") == 0)
        method = su::generalized;
    else {
        return unknown_method;
    }

    std::ifstream ifs(tree_filename);
    std::string content = std::string(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
    su::biom table = su::biom(biom_filename);
    su::BPTree tree = su::BPTree(content);
    
    std::unordered_set<std::string> to_keep(table.obs_ids.begin(), table.obs_ids.end());
    su::BPTree tree_sheared = tree.shear(to_keep).collapse();

    std::vector<double*> dm_stripes; 
    std::vector<double*> dm_stripes_total; 
    dm_stripes.resize((table.n_samples + 1) / 2);
    dm_stripes_total.resize((table.n_samples + 1) / 2);
    
    unsigned int chunksize = dm_stripes.size() / nthreads;
    unsigned int start = 0;
    unsigned int end = dm_stripes.size();
    
    std::vector<su::task_parameters> tasks;
    tasks.resize(nthreads);

    std::vector<std::thread> threads(nthreads);
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        tasks[tid].tid = tid;
        tasks[tid].start = start; // stripe start
        tasks[tid].stop = start + chunksize;  // stripe end 
        tasks[tid].n_samples = table.n_samples;
        tasks[tid].g_unifrac_alpha = alpha;
        start = start + chunksize;
    }
    // the last thread gets any trailing bits
    tasks[threads.size() - 1].stop = end;
    
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        if(variance_adjust)
            threads[tid] = std::thread(su::unifrac_vaw, 
                                       std::ref(table),
                                       std::ref(tree_sheared), 
                                       method, 
                                       std::ref(dm_stripes), 
                                       std::ref(dm_stripes_total), 
                                       &tasks[tid]);
        else
            threads[tid] = std::thread(su::unifrac, 
                                       std::ref(table),
                                       std::ref(tree_sheared), 
                                       method, 
                                       std::ref(dm_stripes), 
                                       std::ref(dm_stripes_total), 
                                       &tasks[tid]);
    }

    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        threads[tid].join();
    }

    initialize_mat(result, table); 
    result->condensed_form = su::deconvolute_stripes(dm_stripes, table.n_samples);
    destroy_stripes(dm_stripes, dm_stripes_total, table.n_samples);

    return okay;
}


