/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "tree.hpp"
#include "biom_interface.hpp"
#include "unifrac.hpp"
#include "affinity.hpp"
#include <unordered_map>
#include <cstdlib>
#include <thread>
#include <signal.h>
#include <stdarg.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>

#include "unifrac_internal.hpp"

// We will always have the CPU version
#define SUCMP_NM  su_cpu
#include "unifrac_cmp.hpp"
#undef SUCMP_NM

#ifdef UNIFRAC_ENABLE_ACC
#define SUCMP_NM  su_acc
#include "unifrac_cmp.hpp"
#undef SUCMP_NM
#endif

using namespace su;

std::string su::test_table_ids_are_subset_of_tree(su::biom_interface &table, su::BPTree &tree) {
    std::unordered_set<std::string> tip_names = tree.get_tip_names();
    std::unordered_set<std::string>::const_iterator hit;
    std::string a_missing_name = "";

    for(auto i : table.obs_ids) {
        hit = tip_names.find(i);
        if(hit == tip_names.end()) {
            a_missing_name = i;
            break;
        }
    }

    return a_missing_name;
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


void su::stripes_to_condensed_form(std::vector<double*> &stripes, uint32_t n, double* cf, unsigned int start, unsigned int stop) {
    // n must be >= 2, but that should be enforced upstream as that would imply
    // computing unifrac on a single sample.

    uint64_t comb_N = comb_2(n);
    for(unsigned int stripe = start; stripe < stop; stripe++) {
        // compute the (i, j) position of each element in each stripe
        uint64_t i = 0;
        uint64_t j = stripe + 1;
        for(uint64_t k = 0; k < n; k++, i++, j++) {
            if(j == n) {
                i = 0;
                j = n - (stripe + 1);
            }
            // determine the position in the condensed form vector for a given (i, j)
            // based off of
            // https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
            uint64_t comb_N_minus_i = comb_2(n - i);
            cf[comb_N - comb_N_minus_i + (j - i - 1)] = stripes[stripe][k];
        }
    }
}


// write in a 2D matrix 
// also suitable for writing to disk
template<class TReal>
void su::condensed_form_to_matrix_T(const double*  __restrict__ cf, const uint32_t n, TReal*  __restrict__ buf2d) {
     const uint64_t comb_N = su::comb_2(n);
     for(uint64_t i = 0; i < n; i++) {
        for(uint64_t j = 0; j < n; j++) {
            TReal v;
            if(i < j) { // upper triangle
                const uint64_t comb_N_minus = su::comb_2(n - i);
                v = cf[comb_N - comb_N_minus + (j - i - 1)];
            } else if (i > j) { // lower triangle
                const uint64_t comb_N_minus = su::comb_2(n - j);
                v = cf[comb_N - comb_N_minus + (i - j - 1)];
            } else {
                v = 0.0;
            }
            buf2d[i*n+j] = v;
        }
     }
}


// make sure it is instantiated
template void su::condensed_form_to_matrix_T<double>(const double*  __restrict__ cf, const uint32_t n, double*  __restrict__ buf2d);
template void su::condensed_form_to_matrix_T<float>(const double*  __restrict__ cf, const uint32_t n, float*  __restrict__ buf2d);

void su::condensed_form_to_matrix(const double*  __restrict__ cf, const uint32_t n, double*  __restrict__ buf2d) {
  su::condensed_form_to_matrix_T<double>(cf,n,buf2d);
}

void su::condensed_form_to_matrix_fp32(const double*  __restrict__ cf, const uint32_t n, float*  __restrict__ buf2d) {
  su::condensed_form_to_matrix_T<float>(cf,n,buf2d);
}

/*
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
 */


// Helper class
// Will cache pointers and automatically release stripes when all elements are used
class OnceManagedStripes {
   private:
    const uint32_t n_samples;
    const uint32_t n_stripes;
    const ManagedStripes &stripes;
    std::vector<const double *> stripe_ptr;
    std::vector<uint32_t> stripe_accessed;

    const double *get_stripe(const uint32_t stripe) {
      if (stripe_ptr[stripe]==0) stripe_ptr[stripe]=stripes.get_stripe(stripe);
      return stripe_ptr[stripe];
    }

    void release_stripe(const uint32_t stripe) {
       stripes.release_stripe(stripe);
       stripe_ptr[stripe]=0;
    }

   public:
    OnceManagedStripes(const ManagedStripes &_stripes, const uint32_t _n_samples, const uint32_t _n_stripes)
    : n_samples(_n_samples), n_stripes(_n_stripes)
    , stripes(_stripes)
    , stripe_ptr(n_stripes)
    , stripe_accessed(n_stripes)
    {}
    
    ~OnceManagedStripes()
    {
      for(uint32_t i = 0; i < n_stripes; i++) {
        if (stripe_ptr[i]!=0) {
           release_stripe(i); 
        }
      }
    }

    double get_val(const uint32_t stripe, const uint32_t el)
    {
      if (stripe_ptr[stripe]==0) stripe_ptr[stripe]=stripes.get_stripe(stripe); 
      const double *mystripe = stripe_ptr[stripe];
      double val = mystripe[el];

      stripe_accessed[stripe]++;
      if (stripe_accessed[stripe]==n_samples) release_stripe(stripe); // we will not use this stripe anymore

      return val;
    }


};

// write in a 2D matrix 
// also suitable for writing to disk
template<class TReal>
void su::stripes_to_matrix_T(const ManagedStripes &_stripes, const uint32_t n_samples, const uint32_t n_stripes, TReal*  __restrict__ buf2d, uint32_t tile_size) {
    // n_samples must be >= 2, but that should be enforced upstream as that would imply
    // computing unifrac on a single sample.

    // tile for  for better memory access pattern
    const uint32_t TILE = (tile_size>0) ? tile_size : (128/sizeof(TReal));
    const uint32_t n_samples_tup = (n_samples+(TILE-1))/TILE; // round up

    OnceManagedStripes stripes(_stripes, n_samples, n_stripes);

    
    for(uint32_t oi = 0; oi < n_samples_tup; oi++) { // off diagonal
      // alternate between inner and outer off-diagonal, due to wrap around in stripes
      const uint32_t o = ((oi%2)==0) ? \
                           (oi/2)*TILE :                  /* close to diagonal */ \
                           (n_samples_tup-(oi/2)-1)*TILE; /* far from diagonal */
 
      for(uint32_t d = 0; d < (n_samples-o); d+=TILE) { // diagonal

         uint32_t iOut = d;
         uint32_t jOut = d+o;

         uint32_t iMax = std::min(iOut+TILE,n_samples);
         uint32_t jMax = std::min(jOut+TILE,n_samples);


        if (iOut==jOut) { 
          // on diagonal
          for(uint64_t i = iOut; i < iMax; i++) {
             buf2d[i*n_samples+i] = 0.0;

             int64_t stripe=0;

             uint64_t j = i+1;
             for(; (stripe<n_stripes) && (j<jMax); stripe++, j++) {
               TReal val = stripes.get_val(stripe, i);
               buf2d[i*n_samples+j] = val;
             }

             if (j<n_samples) { // implies strip==n_stripes, we are really looking at the mirror
               stripe=n_samples-n_stripes-1;
               for(; j < jMax; j++) {
                 --stripe;
                 TReal val = stripes.get_val(stripe, j);
                 buf2d[i*n_samples+j] = val;
               }
             }
          }

          // lower triangle
          for(uint64_t i = iOut+1; i < iMax; i++) {
            for(uint64_t j = jOut; j < i; j++) {
              buf2d[i*n_samples+j] = buf2d[j*n_samples+i];
            }
          }

        } else if (iOut<jOut) {
          // off diagonal
          for(uint64_t i = iOut; i < iMax; i++) {
             unsigned int stripe=0;

             uint64_t j = i+1;
             // we are off diagonal, so adjust
             stripe += (jOut-j);
             j=jOut;
             if (stripe>n_stripes) {
               // ops, we overshoot... roll back
               j-=(stripe-n_stripes);
               stripe=n_stripes;
             }
             for(; (stripe<n_stripes) && (j<jMax); stripe++, j++) {
                TReal val = stripes.get_val(stripe, i);
                buf2d[i*n_samples+j] = val;
             }

             if (j<jMax) { // implies strip==n_stripes, we are really looking at the mirror
               stripe=n_samples-n_stripes-1;
               if (j<jOut) {
                 stripe -= (jOut-j); // note: should not be able to overshoot
                 j=jOut;
               }
               for(; j < jMax; j++) {
                 --stripe;
                 TReal val = stripes.get_val(stripe, j);
                 buf2d[i*n_samples+j] = val;
               }
             }
          }

          // do the other off-diagonal immediately, so it is still in cache
          for(uint64_t j = jOut; j < jMax; j++) {
            for(uint64_t i = iOut; i < iMax; i++) {
              buf2d[j*n_samples+i] = buf2d[i*n_samples+j];
            }
          }
        }
 
      } //for jOut 
    } // for iOut
}

// Make sure it gets instantiated
template void su::stripes_to_matrix_T<double>(const ManagedStripes &stripes, const uint32_t n_samples, const uint32_t n_stripes, double*  __restrict__ buf2d, uint32_t tile_size);
template void su::stripes_to_matrix_T<float>(const ManagedStripes &stripes, const uint32_t n_samples, const uint32_t n_stripes, float*  __restrict__ buf2d, uint32_t tile_size);

void su::stripes_to_matrix(const ManagedStripes &stripes, const uint32_t n_samples, const uint32_t n_stripes, double*  __restrict__ buf2d, uint32_t tile_size) {
   return su::stripes_to_matrix_T<double>(stripes, n_samples, n_stripes, buf2d, tile_size);
}

void su::stripes_to_matrix_fp32(const ManagedStripes &stripes, const uint32_t n_samples, const uint32_t n_stripes, float*  __restrict__ buf2d, uint32_t tile_size) {
  return su::stripes_to_matrix_T<float>(stripes, n_samples, n_stripes, buf2d, tile_size);
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

// Computes Faith's PD for the samples in  `table` over the phylogenetic
// tree given by `tree`.
// Assure that tree does not contain ids that are not in table
void su::faith_pd(biom_interface &table,
                  BPTree &tree,
                  double* result) {
    PropStack<double> propstack(table.n_samples);

    uint32_t node;
    double *node_proportions;
    double length;

    // for node in postorderselect
    for(unsigned int k = 0; k < (tree.nparens / 2) - 1; k++) {
        node = tree.postorderselect(k);
        // get branch length
        length = tree.lengths[node];

        // get node proportions and set intermediate scores
        node_proportions = propstack.pop(node);
        set_proportions(node_proportions, tree, node, table, propstack);

        for (unsigned int sample = 0; sample < table.n_samples; sample++){
            // calculate contribution of node to score
            result[sample] += (node_proportions[sample] > 0) * length;
        }
    }
}


#ifdef UNIFRAC_ENABLE_ACC

// test only once, then use persistent value
static int proc_use_acc = -1;

inline bool use_acc() {
 if (proc_use_acc!=-1) return (proc_use_acc!=0);
 int has_nvidia_gpu_rc = access("/proc/driver/nvidia/gpus", F_OK);

 bool print_info = false;

 if (const char* env_p = std::getenv("UNIFRAC_GPU_INFO")) {
   print_info = true;
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     print_info = false;
   }
 }


 if (has_nvidia_gpu_rc != 0) {
   if (print_info) printf("INFO (unifrac): GPU not found, using CPU\n");
   proc_use_acc=0;
   return false;
 }

 if (const char* env_p = std::getenv("UNIFRAC_USE_GPU")) {
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     if (print_info) printf("INFO (unifrac): Use of GPU explicitly disabled, using CPU\n");
     proc_use_acc=0;
     return false;
   }
 }

 if (print_info) printf("INFO (unifrac): Using GPU\n");
 proc_use_acc=1;
 return true;
}
#endif

void su::unifrac(biom_interface &table,
                 BPTree &tree,
                 Method unifrac_method,
                 std::vector<double*> &dm_stripes,
                 std::vector<double*> &dm_stripes_total,
                 const su::task_parameters* task_p) {
#ifdef UNIFRAC_ENABLE_ACC
  if (use_acc()) {
    su_acc::unifrac(table, tree, unifrac_method, dm_stripes, dm_stripes_total, task_p);
  } else {
#else
  if (true) {
#endif
    su_cpu::unifrac(table, tree, unifrac_method, dm_stripes, dm_stripes_total, task_p);
  }
}


void su::unifrac_vaw(biom_interface &table,
                     BPTree &tree,
                     Method unifrac_method,
                     std::vector<double*> &dm_stripes,
                     std::vector<double*> &dm_stripes_total,
                     const su::task_parameters* task_p) {
#ifdef UNIFRAC_ENABLE_ACC
  if (use_acc()) {
   su_acc::unifrac_vaw(table, tree, unifrac_method, dm_stripes, dm_stripes_total, task_p);
  } else {
#else
  if (true) {
#endif
   su_cpu::unifrac_vaw(table, tree, unifrac_method, dm_stripes, dm_stripes_total, task_p);
  }
}


void su::process_stripes(biom_interface &table,
                         BPTree &tree_sheared,
                         Method method,
                         bool variance_adjust,
                         std::vector<double*> &dm_stripes,
                         std::vector<double*> &dm_stripes_total,
                         std::vector<std::thread> &threads,
                         std::vector<su::task_parameters> &tasks) {

    // register a signal handler so we can ask the master thread for its
    // progress
    register_report_status();

    // cannot use threading with openacc or openmp
    for(unsigned int tid = 0; tid < threads.size(); tid++) {
        if(variance_adjust)
            su::unifrac_vaw(
                                       std::ref(table),
                                       std::ref(tree_sheared),
                                       method,
                                       std::ref(dm_stripes),
                                       std::ref(dm_stripes_total),
                                       &tasks[tid]);
        else
            su::unifrac(
                                       std::ref(table),
                                       std::ref(tree_sheared),
                                       method,
                                       std::ref(dm_stripes),
                                       std::ref(dm_stripes_total),
                                       &tasks[tid]);
    }

    remove_report_status();
}
