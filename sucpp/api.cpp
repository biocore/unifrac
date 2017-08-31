#include "api.hpp"
#include "biom.hpp"
#include "tree.hpp"
#include "unifrac.hpp"
#include <fstream>
#include <iomanip>
#include <thread>
#include <cstring>

#define CHECK_FILE(filename, err) if(!is_file_exists(filename)) { \
                                      return err;                 \
                                  }

#define SET_METHOD(requested_method, err) Method method;                                                       \
                                          if(std::strcmp(requested_method, "unweighted") == 0)                 \
                                              method = unweighted;                                             \
                                          else if(std::strcmp(requested_method, "weighted_normalized") == 0)   \
                                              method = weighted_normalized;                                    \
                                          else if(std::strcmp(requested_method, "weighted_unnormalized") == 0) \
                                              method = weighted_unnormalized;                                  \
                                          else if(std::strcmp(requested_method, "generalized") == 0)           \
                                              method = generalized;                                            \
                                          else {                                                               \
                                              return err;                                                      \
                                          }

using namespace su;
using namespace std;

// https://stackoverflow.com/a/19841704/19741
bool is_file_exists(const char *fileName) {
    std::ifstream infile(fileName);
        return infile.good();
}


void destroy_stripes(vector<double*> &dm_stripes, vector<double*> &dm_stripes_total, unsigned int n_samples, unsigned int stripe_start, unsigned int stripe_stop) {
    unsigned int n_rotations = (n_samples + 1) / 2;

    for(unsigned int i = 0; i < n_rotations; i++) {
        // if a stripe_stop is specified, and if we're in the stripe window, do not free
        // dm_stripes  
        if(stripe_stop == 0)
            free(dm_stripes[i]);
        else if(i < stripe_start && i >= stripe_stop)
            free(dm_stripes[i]);

        if(dm_stripes_total[i] != NULL) {
            free(dm_stripes_total[i]);
        }
    }
}


void initialize_mat(mat_t* &result, biom &table) {
    result = (mat_t*)malloc(sizeof(mat));
    result->n_samples = table.n_samples;
    
    result->cf_size = su::comb_2(table.n_samples); 
    result->is_square = true;  // future support for dissimilarity matrices
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);
    result->condensed_form = (double*)malloc(sizeof(double) * su::comb_2(table.n_samples));

    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table.sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table.sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
    }
}


void initialize_partial_mat(partial_mat_t* &result, biom &table, std::vector<double*> &dm_stripes, unsigned int stripe_start, unsigned int stripe_stop) {
    result = (partial_mat_t*)malloc(sizeof(partial_mat));
    result->n_samples = table.n_samples;
    
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);
    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table.sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table.sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
    }

    result->stripes = (double**)malloc(sizeof(double*) * (stripe_stop - stripe_start));
    result->stripe_start = stripe_start;
    result->stripe_stop = stripe_stop;
    for(unsigned int i = stripe_start; i < stripe_stop; i++) {
        result->stripes[i - stripe_start] = dm_stripes[i];
    }
}


void destroy_mat(mat_t** result) {
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        free((*result)->sample_ids[i]);
    };
    free((*result)->sample_ids);
    free((*result)->condensed_form);
    free(*result);
}

void set_tasks(std::vector<su::task_parameters> &tasks,
               double alpha,
               unsigned int n_samples,
               unsigned int stripe_start, 
               unsigned int stripe_stop, 
               unsigned int nthreads) {

    // compute from start to the max possible stripe if stop doesn't make sense
    if(stripe_stop <= stripe_start)
        stripe_stop = (n_samples + 1) / 2;

    /* chunking strategy is to balance as much as possible. eg if there are 15 stripes
     * and 4 threads, our goal is to assign 4 stripes to 3 threads, and 3 stripes to one thread.
     *
     * we use the remaining the chunksize for bins which cannot be full maximally
     */
    unsigned int fullchunk = ((stripe_stop - stripe_start) + nthreads - 1) / nthreads;  // this computes the ceiling
    unsigned int smallchunk = (stripe_stop - stripe_start) / nthreads;

    unsigned int n_fullbins = (stripe_start - stripe_stop) % nthreads;
    if(n_fullbins == 0)
        n_fullbins = nthreads;

    unsigned int start = stripe_start;
    
    for(unsigned int tid = 0; tid < nthreads; tid++) {
        tasks[tid].tid = tid;
        tasks[tid].start = start; // stripe start

        if(tid < n_fullbins) {
            tasks[tid].stop = start + fullchunk;  // stripe end 
            start = start + fullchunk;
        } else {
            tasks[tid].stop = start + smallchunk;  // stripe end 
            start = start + smallchunk;
        }

        tasks[tid].n_samples = n_samples;
        tasks[tid].g_unifrac_alpha = alpha;
    }
}

compute_status partial(const char* biom_filename, const char* tree_filename, 
                       const char* unifrac_method, bool variance_adjust, double alpha,
                       unsigned int nthreads, unsigned int stripe_start, unsigned int stripe_stop,
                       partial_mat_t** result) {

    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    SET_METHOD(unifrac_method, unknown_method)
    //PARSE_SYNC_TREE_TABLE(tree_filename, table_filename)

    std::vector<su::task_parameters> tasks(nthreads);
    std::vector<std::thread> threads(nthreads);

    //// for testing
    su::biom table = su::biom(biom_filename); 

    // we resize to the largest number of possible stripes even if only computing
    // partial, however we do not allocate arrays for non-computed stripes so
    // there is a little memory waste here but should be on the order of
    // 8 bytes * N samples per vector. 
    std::vector<double*> dm_stripes(stripe_stop); 
    std::vector<double*> dm_stripes_total(stripe_stop);

    set_tasks(tasks, alpha, table.n_samples, stripe_start, stripe_stop, nthreads);
    su::process_stripes(biom_filename, tree_filename, method, variance_adjust, dm_stripes, dm_stripes_total, threads, tasks);

    initialize_partial_mat(*result, table, dm_stripes, stripe_start, stripe_stop);
    destroy_stripes(dm_stripes, dm_stripes_total, table.n_samples, stripe_start, stripe_stop);
    
    return okay;
}


compute_status one_off(const char* biom_filename, const char* tree_filename, 
                       const char* unifrac_method, bool variance_adjust, double alpha,
                       unsigned int nthreads, mat_t** result) {

    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    SET_METHOD(unifrac_method, unknown_method)
    //PARSE_SYNC_TREE_TABLE(tree_filename, table_filename)


    //// for testing
    su::biom table = su::biom(biom_filename); 


    std::vector<su::task_parameters> tasks(nthreads);
    std::vector<std::thread> threads(nthreads);

    // we resize to the largest number of possible stripes even if only computing
    // partial, however we do not allocate arrays for non-computed stripes so
    // there is a little memory waste here but should be on the order of
    // 8 bytes * N samples per vector. 
    std::vector<double*> dm_stripes(comb_2(table.n_samples)); 
    std::vector<double*> dm_stripes_total(comb_2(table.n_samples));

    set_tasks(tasks, alpha, table.n_samples, 0, 0, nthreads);
    su::process_stripes(biom_filename, tree_filename, method, variance_adjust, dm_stripes, dm_stripes_total, threads, tasks);

    initialize_mat(*result, table);
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        threads[tid] = std::thread(su::stripes_to_condensed_form, 
                                   std::ref(dm_stripes), 
                                   table.n_samples,
                                   std::ref((*result)->condensed_form),
                                   tasks[tid].start,
                                   tasks[tid].stop);
    } 
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        threads[tid].join();
    }

    destroy_stripes(dm_stripes, dm_stripes_total, table.n_samples, 0, 0);

    return okay;
}
