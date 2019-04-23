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

#define PARSE_SYNC_TREE_TABLE(tree_filename, table_filename) std::ifstream ifs(tree_filename);                                        \
                                                             std::string content = std::string(std::istreambuf_iterator<char>(ifs),   \
                                                                                               std::istreambuf_iterator<char>());     \
                                                             su::BPTree tree = su::BPTree(content);                                   \
                                                             su::biom table = su::biom(biom_filename);                                \
                                                             std::string bad_id = su::test_table_ids_are_subset_of_tree(table, tree); \
                                                             if(bad_id != "") {                                                       \
                                                                 return table_and_tree_do_not_overlap;                                \
                                                             }                                                                        \
                                                             std::unordered_set<std::string> to_keep(table.obs_ids.begin(),           \
                                                                                                     table.obs_ids.end());            \
                                                             su::BPTree tree_sheared = tree.shear(to_keep).collapse();


using namespace su;
using namespace std;

// https://stackoverflow.com/a/19841704/19741
bool is_file_exists(const char *fileName) {
    std::ifstream infile(fileName);
        return infile.good();
}


void destroy_stripes(vector<double*> &dm_stripes, vector<double*> &dm_stripes_total, unsigned int n_samples, unsigned int stripe_start, unsigned int stripe_stop) {
    unsigned int n_rotations = (n_samples + 1) / 2;

    if(stripe_stop == 0) {
        for(unsigned int i = 0; i < n_rotations; i++) {
            free(dm_stripes[i]);
            if(dm_stripes_total[i] != NULL)
                free(dm_stripes_total[i]);
        }
    } else {
        // if a stripe_stop is specified, and if we're in the stripe window, do not free
        // dm_stripes. this is done as the pointers in dm_stripes are assigned to the partial_mat_t
        // and subsequently freed in destroy_partial_mat. but, we do need to free dm_stripes_total
        // if appropriate
        for(unsigned int i = stripe_start; i < stripe_stop; i++) {
            if(dm_stripes_total[i] != NULL)
                free(dm_stripes_total[i]);
        }
    }
}


void initialize_mat(mat_t* &result, biom &table, bool is_upper_triangle) {
    result = (mat_t*)malloc(sizeof(mat));
    result->n_samples = table.n_samples;

    result->cf_size = su::comb_2(table.n_samples);
    result->is_upper_triangle = is_upper_triangle;
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);
    result->condensed_form = (double*)malloc(sizeof(double) * su::comb_2(table.n_samples));

    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table.sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table.sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
    }
}

void initialize_results_vec(r_vec* &result, biom& table){
    result = (r_vec*)malloc(sizeof(results_vec));
    result->n_samples = table.n_samples;
    result->values = (double*)malloc(sizeof(double) * result->n_samples);
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);

    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table.sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table.sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
        result->values[i] = 0;
    }

}

void initialize_mat_no_biom(mat_t* &result, char** sample_ids, unsigned int n_samples, bool is_upper_triangle) {
    result = (mat_t*)malloc(sizeof(mat));
    result->n_samples = n_samples;

    result->cf_size = su::comb_2(n_samples);
    result->is_upper_triangle = is_upper_triangle;
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);
    result->condensed_form = (double*)malloc(sizeof(double) * su::comb_2(n_samples));

    for(unsigned int i = 0; i < n_samples; i++) {
        result->sample_ids[i] = strdup(sample_ids[i]);
    }
}

void initialize_partial_mat(partial_mat_t* &result, biom &table, std::vector<double*> &dm_stripes, unsigned int stripe_start, unsigned int stripe_stop, bool is_upper_triangle) {
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
    result->is_upper_triangle = is_upper_triangle;
    result->stripe_total = dm_stripes.size();

    for(unsigned int i = stripe_start; i < stripe_stop; i++) {
        result->stripes[i - stripe_start] = dm_stripes[i];
    }
}

void destroy_results_vec(r_vec** result) {
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        free((*result)->sample_ids[i]);
    };
    free((*result)->sample_ids);
    free((*result)->values);
    free(*result);

}
void destroy_mat(mat_t** result) {
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        free((*result)->sample_ids[i]);
    };
    free((*result)->sample_ids);
    free((*result)->condensed_form);
    free(*result);
}

void destroy_partial_mat(partial_mat_t** result) {
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        if((*result)->sample_ids[i] != NULL)
            free((*result)->sample_ids[i]);
    };
    if((*result)->sample_ids != NULL)
        free((*result)->sample_ids);

    unsigned int n_stripes = (*result)->stripe_stop - (*result)->stripe_start;
    for(unsigned int i = 0; i < n_stripes; i++)
        if((*result)->stripes[i] != NULL)
            free((*result)->stripes[i]);
    if((*result)->stripes != NULL)
        free((*result)->stripes);

    free(*result);
}

void set_tasks(std::vector<su::task_parameters> &tasks,
               double alpha,
               unsigned int n_samples,
               unsigned int stripe_start,
               unsigned int stripe_stop,
               bool bypass_tips,
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

    unsigned int n_fullbins = (stripe_stop - stripe_start) % nthreads;
    if(n_fullbins == 0)
        n_fullbins = nthreads;

    unsigned int start = stripe_start;

    for(unsigned int tid = 0; tid < nthreads; tid++) {
        tasks[tid].tid = tid;
        tasks[tid].start = start; // stripe start
        tasks[tid].bypass_tips = bypass_tips;

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
                       const char* unifrac_method, bool variance_adjust, double alpha, bool bypass_tips,
                       unsigned int nthreads, unsigned int stripe_start, unsigned int stripe_stop,
                       partial_mat_t** result) {

    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    SET_METHOD(unifrac_method, unknown_method)
    PARSE_SYNC_TREE_TABLE(tree_filename, table_filename)

    // we resize to the largest number of possible stripes even if only computing
    // partial, however we do not allocate arrays for non-computed stripes so
    // there is a little memory waste here but should be on the order of
    // 8 bytes * N samples per vector.
    std::vector<double*> dm_stripes((table.n_samples + 1) / 2);
    std::vector<double*> dm_stripes_total((table.n_samples + 1) / 2);

    if(nthreads > dm_stripes.size()) {
        fprintf(stderr, "More threads were requested than stripes. Using %d threads.\n");
        nthreads = dm_stripes.size();
    }

    std::vector<su::task_parameters> tasks(nthreads);
    std::vector<std::thread> threads(nthreads);

    if(((table.n_samples + 1) / 2) < stripe_stop) {
        fprintf(stderr, "Stopping stripe is out-of-bounds, max %d\n", (table.n_samples + 1) / 2);
        exit(EXIT_FAILURE);
    }

    set_tasks(tasks, alpha, table.n_samples, stripe_start, stripe_stop, bypass_tips, nthreads);
    su::process_stripes(table, tree_sheared, method, variance_adjust, dm_stripes, dm_stripes_total, threads, tasks);

    initialize_partial_mat(*result, table, dm_stripes, stripe_start, stripe_stop, true);  // true -> is_upper_triangle
    destroy_stripes(dm_stripes, dm_stripes_total, table.n_samples, stripe_start, stripe_stop);

    return okay;
}

compute_status faith_pd_one_off(const char* biom_filename, const char* tree_filename,
                                r_vec** result){
    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    PARSE_SYNC_TREE_TABLE(tree_filename, table_filename)

    initialize_results_vec(*result, table);

    // compute faithpd
    su::faith_pd(table, tree_sheared, std::ref((*result)->values));

    return okay;

}


compute_status one_off(const char* biom_filename, const char* tree_filename,
                       const char* unifrac_method, bool variance_adjust, double alpha,
                       bool bypass_tips, unsigned int nthreads, mat_t** result) {

    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    SET_METHOD(unifrac_method, unknown_method)
    PARSE_SYNC_TREE_TABLE(tree_filename, table_filename)

    // we resize to the largest number of possible stripes even if only computing
    // partial, however we do not allocate arrays for non-computed stripes so
    // there is a little memory waste here but should be on the order of
    // 8 bytes * N samples per vector.
    std::vector<double*> dm_stripes((table.n_samples + 1) / 2);
    std::vector<double*> dm_stripes_total((table.n_samples + 1) / 2);

    if(nthreads > dm_stripes.size()) {
        fprintf(stderr, "More threads were requested than stripes. Using %d threads.\n");
        nthreads = dm_stripes.size();
    }

    std::vector<su::task_parameters> tasks(nthreads);
    std::vector<std::thread> threads(nthreads);

    set_tasks(tasks, alpha, table.n_samples, 0, 0, bypass_tips, nthreads);
    su::process_stripes(table, tree_sheared, method, variance_adjust, dm_stripes, dm_stripes_total, threads, tasks);

    initialize_mat(*result, table, true);  // true -> is_upper_triangle
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

IOStatus write_mat(const char* output_filename, mat_t* result) {
    std::ofstream output;
    output.open(output_filename);

    uint64_t comb_N = su::comb_2(result->n_samples);
    uint64_t comb_N_minus = 0;
    double v;

    for(unsigned int i = 0; i < result->n_samples; i++)
        output << "\t" << result->sample_ids[i];
    output << std::endl;

    for(unsigned int i = 0; i < result->n_samples; i++) {
        output << result->sample_ids[i];
        for(unsigned int j = 0; j < result->n_samples; j++) {
            if(i < j) { // upper triangle
                comb_N_minus = su::comb_2(result->n_samples - i);
                v = result->condensed_form[comb_N - comb_N_minus + (j - i - 1)];
            } else if (i > j) { // lower triangle
                comb_N_minus = su::comb_2(result->n_samples - j);
                v = result->condensed_form[comb_N - comb_N_minus + (i - j - 1)];
            } else {
                v = 0.0;
            }
            output << std::setprecision(16) << "\t" << v;
        }
        output << std::endl;
    }

    return write_okay;
}

IOStatus write_partial(const char* output_filename, partial_mat_t* result) {
    std::ofstream output;
    output.open(output_filename, std::ios::binary);
    if(!output.is_open())
        return open_error;

    uint32_t n_stripes = result->stripe_stop - result->stripe_start;
    std::string magic(PARTIAL_MAGIC);
    uint32_t magic_len = magic.length();

    /* header information */
    output.write(reinterpret_cast<const char*>(&magic_len),                 sizeof(uint16_t));
    output << magic;
    output.write(reinterpret_cast<const char*>(&result->n_samples),         sizeof(uint32_t));
    output.write(reinterpret_cast<const char*>(&n_stripes),                 sizeof(uint32_t));
    output.write(reinterpret_cast<const char*>(&result->stripe_start),      sizeof(uint32_t));
    output.write(reinterpret_cast<const char*>(&result->stripe_total),      sizeof(uint32_t));
    output.write(reinterpret_cast<const char*>(&result->is_upper_triangle), sizeof(uint8_t));

    /* sample IDs */
    for(unsigned int i = 0; i < result->n_samples; i++) {
        uint16_t length = strlen(result->sample_ids[i]);
        output.write(reinterpret_cast<const char*>(&length), sizeof(uint16_t));
        output << result->sample_ids[i];
    }

    /* stripe information */
    for(unsigned int i = 0; i < n_stripes; i++) {
        /// :( streamsize didn't seem to work. probably a fancy way to do this, but the regular loop is fast too
        //output.write(reinterpret_cast<const char*>(&result->stripes[i]), std::streamsize(sizeof(double) * result->n_samples));
        for(unsigned int j = 0; j < result->n_samples; j++)
            output.write(reinterpret_cast<const char*>(&result->stripes[i][j]), sizeof(double));
    }

    /* footer */
    output << magic;
    output.close();

    return write_okay;
}

IOStatus _is_partial_file(const char* input_filename) {
    std::ifstream input;
    input.open(input_filename, std::ios::in | std::ios::binary);
    if(!input.is_open())
        return open_error;

    char magic[32];
    uint16_t magic_len;

    input.read((char*)&magic_len, 2);

    // if the length of the magic is unexpected then bail
    if(magic_len <= 0 || magic_len > 32) {
        return magic_incompatible;
    }

    input.read(magic, magic_len);
    if(strncmp(magic, PARTIAL_MAGIC, magic_len) != 0) {
        return magic_incompatible;
    }

    input.close();
    return read_okay;
}

IOStatus read_partial(const char* input_filename, partial_mat_t** result_out) {
    IOStatus err = _is_partial_file(input_filename);

    if(err != read_okay)
        return err;

    std::ifstream input;
    input.open(input_filename, std::ios::binary);

    /* load header */
    uint16_t magic_len;
    input.read((char*)&magic_len, 2);  // magic length

    char header_magic[32];
    input.read(header_magic, magic_len);  // magic
    header_magic[magic_len] = '\0';

    uint32_t n_samples;
    input.read((char*)&n_samples, 4);  // number of samples

    uint32_t n_stripes;
    input.read((char*)&n_stripes, 4);  // number of stripes

    uint32_t stripe_start;
    input.read((char*)&stripe_start, 4);  // stripe start

    uint32_t stripe_total;
    input.read((char*)&stripe_total, 4);  // stripe total

    bool is_upper_triangle;
    input.read((char*)&is_upper_triangle, 1);  // is_upper_triangle

    /* sanity check header */
    if(n_samples <= 0 || n_stripes <= 0 || stripe_start < 0 || stripe_total <= 0 || is_upper_triangle < 0)
        return bad_header;
    if(stripe_total >= n_samples || n_stripes > stripe_total || stripe_start >= stripe_total || stripe_start + n_stripes > stripe_total)
        return bad_header;

    /* initialize the partial result structure */
    partial_mat_t* result = (partial_mat_t*)malloc(sizeof(partial_mat));
    result->n_samples = n_samples;
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);
    result->stripes = (double**)malloc(sizeof(double*) * (n_stripes));
    result->stripe_start = stripe_start;
    result->stripe_stop = stripe_start + n_stripes;
    result->is_upper_triangle = is_upper_triangle;
    result->stripe_total = stripe_total;

    /* load samples */
    for(int i = 0; i < n_samples; i++) {
        uint16_t sample_length;
        input.read((char*)&sample_length, 2);
        result->sample_ids[i] = (char*)malloc(sizeof(char) * (sample_length + 1));
        input.read(result->sample_ids[i], sample_length);
        result->sample_ids[i][sample_length] = '\0';
    }

    /* load stripes */
    int current_to_load;
    void *ptr;
    for(int i = 0; i < n_stripes; i++) {
        ptr = malloc(sizeof(double) * n_samples);
        if(ptr == NULL) {
            fprintf(stderr, "failed\n");
            exit(1);
        }
        result->stripes[i] = (double*)ptr;
        input.read(reinterpret_cast<char*>(result->stripes[i]), sizeof(double) * n_samples);
    }

    /* sanity check the footer */
    char footer_magic[32];
    input.read(footer_magic, magic_len);
    footer_magic[magic_len] = '\0';

    if(strcmp(header_magic, footer_magic) != 0) {
        return magic_incompatible;
    }

    (*result_out) = result;
    return read_okay;
}

MergeStatus merge_partial(partial_mat_t** partial_mats, int n_partials, unsigned int nthreads, mat_t** result) {
    if(n_partials <= 0) {
        fprintf(stderr, "Zero or less partials.\n");
        exit(EXIT_FAILURE);
    }

    // sanity check
    int n_samples = partial_mats[0]->n_samples;
    bool *stripe_map = (bool*)calloc(sizeof(bool), partial_mats[0]->stripe_total);
    int stripe_count = 0;

    for(int i = 0; i < n_partials; i++) {
        if(partial_mats[i]->n_samples != n_samples) {
            free(stripe_map);
            return partials_mismatch;
        }

        if(partial_mats[0]->stripe_total != partial_mats[i]->stripe_total) {
            free(stripe_map);
            return partials_mismatch;
        }
        if(partial_mats[0]->is_upper_triangle != partial_mats[i]->is_upper_triangle) {
            free(stripe_map);
            return square_mismatch;
        }
        for(int j = 0; j < n_samples; j++) {
            if(strcmp(partial_mats[0]->sample_ids[j], partial_mats[i]->sample_ids[j]) != 0) {
                free(stripe_map);
                return sample_id_consistency;
            }
        }
        for(int j = partial_mats[i]->stripe_start; j < partial_mats[i]->stripe_stop; j++) {
            if(stripe_map[j]) {
                free(stripe_map);
                return stripes_overlap;
            }
            stripe_map[j] = true;
            stripe_count += 1;
        }
    }
    free(stripe_map);

    if(stripe_count != partial_mats[0]->stripe_total) {
        return incomplete_stripe_set;
    }

    std::vector<double*> stripes(partial_mats[0]->stripe_total);
    std::vector<double*> stripes_totals(partial_mats[0]->stripe_total);  // not actually used but destroy_stripes needs this to "exist"
    for(int i = 0; i < n_partials; i++) {
        int n_stripes = partial_mats[i]->stripe_stop - partial_mats[i]->stripe_start;
        for(int j = 0; j < n_stripes; j++) {
            // as this is potentially a large amount of memory, don't copy, just adopt
            *&(stripes[j + partial_mats[i]->stripe_start]) = partial_mats[i]->stripes[j];
        }
    }

    if(nthreads > stripes.size()) {
        fprintf(stderr, "More threads were requested than stripes. Using %d threads.\n");
        nthreads = stripes.size();
    }

    std::vector<su::task_parameters> tasks(nthreads);
    std::vector<std::thread> threads(nthreads);

    initialize_mat_no_biom(*result, partial_mats[0]->sample_ids, n_samples, partial_mats[0]->is_upper_triangle);
    su::stripes_to_condensed_form(stripes, n_samples, (*result)->condensed_form, 0, partial_mats[0]->stripe_total);

    destroy_stripes(stripes, stripes_totals, n_samples, 0, n_partials);

    return merge_okay;
}
