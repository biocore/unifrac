#include "api.hpp"
#include "biom.hpp"
#include "tree.hpp"
#include "unifrac.hpp"
#include <fstream>
#include <iomanip>
#include <thread>
#include <cstring>

#include <fcntl.h>
#include <unistd.h>
#include <lz4.h>


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
                                          else if(std::strcmp(requested_method, "unweighted_fp32") == 0)            \
                                              method = unweighted_fp32;                                             \
                                          else if(std::strcmp(requested_method, "weighted_normalized_fp32") == 0)   \
                                              method = weighted_normalized_fp32;                                    \
                                          else if(std::strcmp(requested_method, "weighted_unnormalized_fp32") == 0) \
                                              method = weighted_unnormalized_fp32;                                  \
                                          else if(std::strcmp(requested_method, "generalized_fp32") == 0)           \
                                              method = generalized_fp32;                                            \
                                          else {                                                               \
                                              return err;                                                      \
                                          }

#define PARSE_SYNC_TREE_TABLE(tree_filename, table_filename) std::ifstream ifs(tree_filename);                                        \
                                                             std::string content = std::string(std::istreambuf_iterator<char>(ifs),   \
                                                                                               std::istreambuf_iterator<char>());     \
                                                             su::BPTree tree = su::BPTree(content);                                   \
                                                             su::biom table = su::biom(biom_filename);                                \
                                                             if(table.n_samples <= 0 | table.n_obs <= 0) {                            \
                                                                 return table_empty;                                                  \
                                                             }                                                                        \
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


void destroy_stripes(vector<double*> &dm_stripes, vector<double*> &dm_stripes_total, unsigned int n_samples,
                     unsigned int stripe_start, unsigned int stripe_stop) {
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
    // Stores results for Faith PD
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

void initialize_partial_mat(partial_mat_t* &result, biom &table, std::vector<double*> &dm_stripes,
                            unsigned int stripe_start, unsigned int stripe_stop, bool is_upper_triangle) {
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
    // for Faith PD
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
        fprintf(stderr, "More threads were requested than stripes. Using %d threads.\n", dm_stripes.size());
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
        fprintf(stderr, "More threads were requested than stripes. Using %d threads.\n", dm_stripes.size());
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
    output.close();

    return write_okay;
}

herr_t write_hdf5_string(hid_t output_file_id,const char *dname, const char *str)
{
  // this is the convoluted way to store a string
  // Will use the FORTRAN forma, so we do not depend on null termination
  hid_t filetype_id = H5Tcopy (H5T_FORTRAN_S1);
  H5Tset_size(filetype_id, strlen(str));
  hid_t memtype_id = H5Tcopy (H5T_C_S1);
  H5Tset_size(memtype_id, strlen(str)+1);

  hsize_t  dims[1] = {1};
  hid_t dataspace_id = H5Screate_simple (1, dims, NULL);

  hid_t dataset_id = H5Dcreate(output_file_id,dname, filetype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT);
  herr_t status = H5Dwrite(dataset_id, memtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, str);

  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Tclose(memtype_id);
  H5Tclose(filetype_id);

  return status;
}

// Internal: Make sure TReal and real_id match
template<class TReal>
IOStatus write_mat_hdf5_D(const char* output_filename, mat_t* result,hid_t real_id, unsigned int compress_level) {
   /* Create a new file using default properties. */
   hid_t output_file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   if (output_file_id<0) return open_error;

   // simple header
   if (write_hdf5_string(output_file_id,"format","BDSM")<0) {
       H5Fclose (output_file_id);
       return open_error;
   }
   if (write_hdf5_string(output_file_id,"version","2020.06")<0) {
       H5Fclose (output_file_id);
       return open_error;
   }

   // save the ids
   {
     hsize_t     dims[1];
     dims[0] = result->n_samples;
     hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

     // this is the convoluted way to store an array of strings
     hid_t datatype_id = H5Tcopy(H5T_C_S1);
     H5Tset_size(datatype_id,H5T_VARIABLE);

     hid_t dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
     if (H5Pset_deflate(dcpl_id, compress_level)<0) return open_error; // just abort on error

     hsize_t     chunks[1];
     chunks[0] = result->n_samples;

     if (H5Pset_chunk (dcpl_id, 1, chunks)) return open_error; // just abort on error

     hid_t dataset_id = H5Dcreate1(output_file_id, "order", datatype_id, dataspace_id, dcpl_id);

     herr_t status = H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, result->sample_ids);

     H5Dclose(dataset_id);
     H5Tclose(datatype_id);
     H5Sclose(dataspace_id);
     H5Pclose(dcpl_id);

     // check status after cleanup, for simplicity
     if (status<0) {
       H5Fclose (output_file_id);
       return open_error;
     }
   }

   // save the matrix
   {
     const uint64_t n_samples = result->n_samples;
     TReal *buf2d = (TReal*) malloc(n_samples*n_samples*sizeof(TReal));
     if (buf2d==NULL) {
       H5Fclose (output_file_id);
       return open_error; // we don't have a better error code
     }

     // but first compute the values to save
     {
       const uint64_t comb_N = su::comb_2(n_samples);
       for(uint64_t i = 0; i < n_samples; i++) {
        for(uint64_t j = 0; j < n_samples; j++) {
            TReal v;
            if(i < j) { // upper triangle
                const uint64_t comb_N_minus = su::comb_2(n_samples - i);
                v = result->condensed_form[comb_N - comb_N_minus + (j - i - 1)];
            } else if (i > j) { // lower triangle
                const uint64_t comb_N_minus = su::comb_2(n_samples - j);
                v = result->condensed_form[comb_N - comb_N_minus + (i - j - 1)];
            } else {
                v = 0.0;
            }
            buf2d[i*n_samples+j] = v;
        }
       }
     }

     hsize_t     dims[2];
     dims[0] = result->n_samples;
     dims[1] = result->n_samples;
     hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

     hid_t dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
     if (H5Pset_deflate(dcpl_id, compress_level)<0) return open_error; // just abort on error

     // shoot for a 0.75M chunk size at double, to fit in default cache
     hsize_t     chunks[2];
     chunks[0] = 1;
     chunks[1] = 96*1024;
     if ( chunks[1]>(result->n_samples) ) {
       chunks[1] = result->n_samples;
       chunks[0] = 96*1024/chunks[1];
       if ( chunks[0]>(result->n_samples) ) {
          chunks[0] = result->n_samples;
       }
     }

     if (H5Pset_chunk (dcpl_id, 2, chunks)) return open_error; // just abort on error

     hid_t dataset_id = H5Dcreate2(output_file_id, "matrix",real_id, dataspace_id,
                                   H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
     herr_t status = H5Dwrite(dataset_id, real_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              buf2d);

     H5Pclose(dcpl_id);
     H5Dclose(dataset_id);
     H5Sclose(dataspace_id);
     free(buf2d);

     // check status after cleanup, for simplicity
     if (status<0) {
       H5Fclose (output_file_id);
       return open_error;
     }
   } 

   H5Fclose (output_file_id);
   return write_okay;
}

IOStatus write_mat_hdf5(const char* output_filename, mat_t* result) {
  return write_mat_hdf5_D<double>(output_filename,result,H5T_IEEE_F64LE,0);
}

IOStatus write_mat_hdf5_fp32(const char* output_filename, mat_t* result) {
  return write_mat_hdf5_D<float>(output_filename,result,H5T_IEEE_F32LE,0);
}

IOStatus write_mat_hdf5_compressed(const char* output_filename, mat_t* result, unsigned int compress_level) {
  return write_mat_hdf5_D<double>(output_filename,result,H5T_IEEE_F64LE,compress_level);
}

IOStatus write_mat_hdf5_fp32_compressed(const char* output_filename, mat_t* result, unsigned int compress_level) {
  return write_mat_hdf5_D<float>(output_filename,result,H5T_IEEE_F32LE,compress_level);
}

IOStatus write_vec(const char* output_filename, r_vec* result) {
    std::ofstream output;
    output.open(output_filename);

    // write sample ids in first column of file and faith's pd in second column
    output << "#SampleID\tfaith_pd" << std::endl;
    for(unsigned int i = 0; i < result->n_samples; i++) {
        output << result->sample_ids[i];
        output << std::setprecision(16) << "\t" << result->values[i];
        output << std::endl;
    }
    output.close();

    return write_okay;
}

IOStatus write_partial(const char* output_filename, partial_mat_t* result) {
    int fd = open(output_filename, O_WRONLY | O_CREAT | O_TRUNC,  S_IRUSR |  S_IWUSR );
    if (fd==-1) return open_error;

    int cnt = -1;

    uint32_t n_stripes = result->stripe_stop - result->stripe_start;

    uint32_t sample_id_length = 0;
    for(unsigned int i = 0; i < result->n_samples; i++) {
        sample_id_length += strlen(result->sample_ids[i])+1;
    }

    {
      char * const samples_buf = (char *)malloc(sample_id_length);
 
      char *samples_ptr = samples_buf;

      /* sample IDs */
      for(unsigned int i = 0; i < result->n_samples; i++) {
          uint32_t length = strlen(result->sample_ids[i])+1;
          memcpy(samples_ptr,result->sample_ids[i],length);
          samples_ptr+= length;
      }

      int max_compressed = LZ4_compressBound(sample_id_length);
      char * const cmp_buf = (char *)malloc(max_compressed);

      int sample_id_length_compressed = LZ4_compress_default(samples_buf,cmp_buf,sample_id_length,max_compressed);
      if (sample_id_length_compressed<1)  {close(fd); return open_error;}

      uint32_t header[8];
      header[0] = PARTIAL_MAGIC_V2;
      header[1] = result->n_samples;
      header[2] = n_stripes;
      header[3] = result->stripe_start;
      header[4] = result->stripe_total;
      header[5] = result->is_upper_triangle;
      header[6] = sample_id_length;
      header[7] = sample_id_length_compressed;

      cnt=write(fd,header, 8 * sizeof(uint32_t));
      if (cnt<1)  {close(fd); return open_error;}

      cnt=write(fd,cmp_buf, sample_id_length_compressed);
      if (cnt<1)  {close(fd); return open_error;}

      free(cmp_buf);
      free(samples_buf);
    }

    {
      int max_compressed = LZ4_compressBound(sizeof(double) * result->n_samples);
      char * const cmp_buf_raw = (char *)malloc(max_compressed+sizeof(uint32_t));
      char * const cmp_buf = cmp_buf_raw + sizeof(uint32_t);

      /* stripe information */
      for(unsigned int i = 0; i < n_stripes; i++) {
        int cmp_size = LZ4_compress_default((const char *) result->stripes[i],cmp_buf,sizeof(double) * result->n_samples,max_compressed);
        if (cmp_size<1)  {close(fd); return open_error;}

        uint32_t *cmp_buf_size_p = (uint32_t *)cmp_buf_raw;
        *cmp_buf_size_p = cmp_size;

        cnt=write(fd, cmp_buf_raw, cmp_size+sizeof(uint32_t));
        if (cnt<1) {return open_error;}
      }

      free(cmp_buf_raw);
    }

    /* footer */
    {
      uint32_t header[1];
      header[0] = PARTIAL_MAGIC_V2;

      cnt=write(fd,header, 1 * sizeof(uint32_t));
      if (cnt<1)  {close(fd); return open_error;}
    }

    close(fd);

    return write_okay;
}

IOStatus _is_partial_file(const char* input_filename) {
    int fd = open(input_filename, O_RDONLY );
    if (fd==-1) return open_error;

    uint32_t header[1];
    int cnt = read(fd,header,sizeof(uint32_t));
    close(fd);

    if (cnt!=sizeof(uint32_t)) return magic_incompatible;
    if ( header[0] != PARTIAL_MAGIC_V2) return magic_incompatible;

    return read_okay;
}

IOStatus read_partial(const char* input_filename, partial_mat_t** result_out) {
    int fd = open(input_filename, O_RDONLY );
    if (fd==-1) return open_error;

    int cnt=-1;

    uint32_t header[8];
    cnt = read(fd,header,8*sizeof(uint32_t));
    if (cnt != (8*sizeof(uint32_t))) {close(fd); return magic_incompatible;}

    if ( header[0] != PARTIAL_MAGIC_V2) {close(fd); return magic_incompatible;}

    const uint32_t n_samples = header[1];
    const uint32_t n_stripes = header[2];
    const uint32_t stripe_start = header[3];
    const uint32_t stripe_total = header[4];
    const bool is_upper_triangle = header[5];

    /* sanity check header */
    if(n_samples <= 0 || n_stripes <= 0 || stripe_total <= 0 || is_upper_triangle < 0)
         {close(fd); return bad_header;}
    if(stripe_total >= n_samples || n_stripes > stripe_total || stripe_start >= stripe_total || stripe_start + n_stripes > stripe_total)
         {close(fd); return bad_header;}

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
    {
      const uint32_t sample_id_length = header[6];
      const uint32_t sample_id_length_compressed = header[7];

      /* sanity check header */
      if (sample_id_length<=0 || sample_id_length_compressed <=0)
         {close(fd); return bad_header;}

      char * const cmp_buf = (char *)malloc(sample_id_length_compressed);
      cnt = read(fd,cmp_buf,sample_id_length_compressed);
      if (cnt != sample_id_length_compressed) {close(fd); return magic_incompatible;}

      char *samples_buf = (char *)malloc(sample_id_length);

      cnt = LZ4_decompress_safe(cmp_buf,samples_buf,sample_id_length_compressed,sample_id_length);
      if (cnt!=sample_id_length) {close(fd); return magic_incompatible;}

      const char *samples_ptr = samples_buf;

      for(int i = 0; i < n_samples; i++) {
        uint32_t sample_length = strlen(samples_ptr);
        if ((samples_ptr+sample_length+1)>(samples_buf+sample_id_length)) {close(fd); return magic_incompatible;}

        result->sample_ids[i] = (char*)malloc(sample_length + 1);
        memcpy(result->sample_ids[i],samples_ptr,sample_length + 1);
        samples_ptr += sample_length + 1;
      }
      free(samples_buf);
      free(cmp_buf);
    }

    /* load stripes */
    {
      int max_compressed = LZ4_compressBound(sizeof(double) * result->n_samples);
      char * const cmp_buf = (char *)malloc(max_compressed+sizeof(uint32_t));

      uint32_t *cmp_buf_size_p = (uint32_t *)cmp_buf;

      cnt = read(fd,cmp_buf_size_p , sizeof(uint32_t) );
      if (cnt != sizeof(uint32_t) ) {close(fd); return magic_incompatible;}

      for(int i = 0; i < n_stripes; i++) {
        uint32_t cmp_size = *cmp_buf_size_p;

        uint32_t read_size = cmp_size;
        if ( (i+1)<n_stripes ) read_size += sizeof(uint32_t); // last one does not have the cmp_size

        cnt = read(fd,cmp_buf , read_size );
        if (cnt != read_size) {close(fd); return magic_incompatible;}

        result->stripes[i] = (double *) malloc(sizeof(double) * n_samples);
        if(result->stripes[i] == NULL) {
            fprintf(stderr, "failed\n");
            exit(1);
        }
        cnt = LZ4_decompress_safe(cmp_buf, (char *) result->stripes[i],cmp_size,sizeof(double) * n_samples);
        if (cnt != ( sizeof(double) * n_samples ) ) {close(fd); return magic_incompatible;}

        cmp_buf_size_p = (uint32_t *)(cmp_buf+cmp_size);
      }

      free(cmp_buf);
    }

    /* sanity check the footer */
    header[0] = 0;
    cnt = read(fd,header,sizeof(uint32_t));
    if (cnt != (sizeof(uint32_t))) {close(fd); return magic_incompatible;}

    if ( header[0] != PARTIAL_MAGIC_V2) {close(fd); return magic_incompatible;}

    close(fd);

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
        fprintf(stderr, "More threads were requested than stripes. Using %d threads.\n", stripes.size());
        nthreads = stripes.size();
    }

    std::vector<su::task_parameters> tasks(nthreads);
    std::vector<std::thread> threads(nthreads);

    initialize_mat_no_biom(*result, partial_mats[0]->sample_ids, n_samples, partial_mats[0]->is_upper_triangle);
    su::stripes_to_condensed_form(stripes, n_samples, (*result)->condensed_form, 0, partial_mats[0]->stripe_total);

    destroy_stripes(stripes, stripes_totals, n_samples, 0, n_partials);

    return merge_okay;
}
