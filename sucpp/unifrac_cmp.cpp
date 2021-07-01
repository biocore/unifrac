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
#include <unordered_map>
#include <cstdlib>
#include <thread>
#include <algorithm>

#include "unifrac_internal.hpp"

#include "unifrac_task.hpp"
// Note: unifrac_task.hpp defines SUCMP_NM, needed by unifrac_cmp.hpp
#include "unifrac_cmp.hpp"

// embed in this file, to properly instantiate the templatized functions
#include "unifrac_task.cpp"

using namespace SUCMP_NM;

template<class TFloat>
inline void initialize_sample_counts(TFloat*& _counts, const su::task_parameters* task_p, const su::biom_interface &table) {
    const unsigned int n_samples = task_p->n_samples;
    const uint64_t  n_samples_r = ((n_samples + UNIFRAC_BLOCK-1)/UNIFRAC_BLOCK)*UNIFRAC_BLOCK; // round up
    TFloat * counts = NULL;
    int err = 0;
    err = posix_memalign((void **)&counts, 4096, sizeof(TFloat) * n_samples_r);
    if(counts == NULL || err != 0) {
        fprintf(stderr, "Failed to allocate %zd bytes, err %d; [%s]:%d\n",
                sizeof(TFloat) * n_samples_r, err, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    for(unsigned int i = 0; i < n_samples; i++) {
        counts[i] = table.sample_counts[i];
    }
   // avoid NaNs
   for(unsigned int i = n_samples; i < n_samples_r; i++) {
       counts[i] = 0.0;
   }

   _counts=counts;
}

template<class TaskT, class TFloat>
inline void unifracTT(const su::biom_interface &table,
                      const su::BPTree &tree,
                      const bool want_total,
                      std::vector<double*> &dm_stripes,
                      std::vector<double*> &dm_stripes_total,
                      const su::task_parameters* task_p) {
    int err;
    // no processor affinity whenusing openacc or openmp

    if(table.n_samples != task_p->n_samples) {
        fprintf(stderr, "Task and table n_samples not equal\n");
        exit(EXIT_FAILURE);
    }
    const unsigned int n_samples = task_p->n_samples;
    const uint64_t  n_samples_r = ((n_samples + UNIFRAC_BLOCK-1)/UNIFRAC_BLOCK)*UNIFRAC_BLOCK; // round up


    su::PropStackMulti<TFloat> propstack_multi(table.n_samples);

    const unsigned int max_emb =  TaskT::RECOMMENDED_MAX_EMBS;

    su::initialize_stripes(std::ref(dm_stripes), std::ref(dm_stripes_total), want_total, task_p);

    TaskT taskObj(std::ref(dm_stripes), std::ref(dm_stripes_total),max_emb,task_p);

    TFloat *lengths = NULL;
    err = posix_memalign((void **)&lengths, 4096, sizeof(TFloat) * max_emb);
    if(err != 0) {
        fprintf(stderr, "posix_memalign(%d) failed: %d\n", sizeof(TFloat) * max_emb,  err);
        exit(EXIT_FAILURE);
    }
#pragma acc enter data create(lengths[:max_emb])

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

    unsigned int k = 0; // index in tree
    const unsigned int max_k = (tree.nparens / 2) - 1;

    const unsigned int num_prop_chunks = propstack_multi.get_num_stacks();
    while (k<max_k) {
          const unsigned int k_start = k;
          unsigned int filled_emb = 0;

          // chunk the progress to maximize cache reuse
#pragma omp parallel for 
          for (unsigned int ck=0; ck<num_prop_chunks; ck++) {
            su::PropStack<TFloat> &propstack = propstack_multi.get_prop_stack(ck);
            const unsigned int tstart = propstack_multi.get_start(ck);
            const unsigned int tend = propstack_multi.get_end(ck);
            unsigned int my_filled_emb = 0;
            unsigned int my_k=k_start;

            while ((my_filled_emb<max_emb) && (my_k<max_k)) {
              const uint32_t node = tree.postorderselect(my_k);
              my_k++;

              TFloat *node_proportions = propstack.pop(node);
              su::set_proportions_range(node_proportions, tree, node, table, tstart, tend, propstack);

              if(task_p->bypass_tips && tree.isleaf(node))
                  continue;

              if (ck==0) { // they all do the same thing, so enough for the first to update the global state
                lengths[filled_emb] = tree.lengths[node];
                filled_emb++;
              }
              taskObj.embed_proportions_range(node_proportions, tstart, tend, my_filled_emb);
              my_filled_emb++;
            }
            if (ck==0) { // they all do the same thing, so enough for the first to update the global state
              k=my_k;
            }
          }

          taskObj.sync_embedded_proportions(filled_emb);
#ifdef _OPENACC
          // lengths may be still in use in async mode, wait
#pragma acc wait
#pragma acc update device(lengths[:filled_emb])
#endif
          taskObj._run(filled_emb,lengths);
          filled_emb=0;

          su::try_report(task_p, k, max_k);
    }

#pragma acc wait

    if(want_total) {
        const uint64_t start_idx = task_p->start;
        const uint64_t stop_idx = task_p->stop;

        TFloat * const dm_stripes_buf = taskObj.dm_stripes.buf;
        const TFloat * const dm_stripes_total_buf = taskObj.dm_stripes_total.buf;

#pragma acc parallel loop collapse(2) present(dm_stripes_buf,dm_stripes_total_buf)
        for(uint64_t i = start_idx; i < stop_idx; i++)
            for(uint64_t j = 0; j < n_samples; j++) {
                uint64_t idx = (i-start_idx)*n_samples_r+j;
                dm_stripes_buf[idx]=dm_stripes_buf[idx]/dm_stripes_total_buf[idx];
                // taskObj.dm_stripes[i][j] = taskObj.dm_stripes[i][j] / taskObj.dm_stripes_total[i][j];
            }
        
    }

#pragma acc exit data delete(lengths[:max_emb])
    free(lengths);
}

void SUCMP_NM::unifrac(const su::biom_interface &table,
                       const su::BPTree &tree,
                       su::Method unifrac_method,
                       std::vector<double*> &dm_stripes,
                       std::vector<double*> &dm_stripes_total,
                        const su::task_parameters* task_p) {
    switch(unifrac_method) {
        case su::unweighted:
            unifracTT<SUCMP_NM::UnifracUnweightedTask<double>,double>(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_normalized:
            unifracTT<SUCMP_NM::UnifracNormalizedWeightedTask<double>,double>(   table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_unnormalized:
            unifracTT<SUCMP_NM::UnifracUnnormalizedWeightedTask<double>,double>( table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::generalized:
            unifracTT<SUCMP_NM::UnifracGeneralizedTask<double>,double>(          table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::unweighted_fp32:
            unifracTT<SUCMP_NM::UnifracUnweightedTask<float >,float>(            table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_normalized_fp32:
            unifracTT<SUCMP_NM::UnifracNormalizedWeightedTask<float >,float>(    table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_unnormalized_fp32:
            unifracTT<SUCMP_NM::UnifracUnnormalizedWeightedTask<float >,float>(  table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::generalized_fp32:
            unifracTT<SUCMP_NM::UnifracGeneralizedTask<float >,float>(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        default:
            fprintf(stderr, "Unknown unifrac task\n");
            exit(1);
            break;
    }
}


template<class TaskT, class TFloat>
inline void unifrac_vawTT(const su::biom_interface &table,
                          const su::BPTree &tree,
                          const bool want_total,
                          std::vector<double*> &dm_stripes,
                          std::vector<double*> &dm_stripes_total,
                          const su::task_parameters* task_p) {
    int err;
    // no processor affinity whenusing openacc or openmp

    if(table.n_samples != task_p->n_samples) {
        fprintf(stderr, "Task and table n_samples not equal\n");
        exit(EXIT_FAILURE);
    }
    const unsigned int n_samples = task_p->n_samples;
    const uint64_t  n_samples_r = ((n_samples + UNIFRAC_BLOCK-1)/UNIFRAC_BLOCK)*UNIFRAC_BLOCK; // round up

    su::PropStackMulti<TFloat> propstack_multi(table.n_samples);
    su::PropStackMulti<TFloat> countstack_multi(table.n_samples);

    const unsigned int max_emb = TaskT::RECOMMENDED_MAX_EMBS;

    TFloat *sample_total_counts;

    initialize_sample_counts(sample_total_counts, task_p, table);
#pragma acc enter data copyin(sample_total_counts[:n_samples_r])
    su::initialize_stripes(std::ref(dm_stripes), std::ref(dm_stripes_total), want_total, task_p);

    TaskT taskObj(std::ref(dm_stripes), std::ref(dm_stripes_total), sample_total_counts, max_emb, task_p);

    TFloat *lengths = NULL;
    err = posix_memalign((void **)&lengths, 4096, sizeof(TFloat) * max_emb);
    if(err != 0) {
        fprintf(stderr, "posix_memalign(%d) failed: %d\n", sizeof(TFloat) * max_emb, err);
        exit(EXIT_FAILURE);
    }
#pragma acc enter data create(lengths[:max_emb])

    unsigned int k = 0; // index in tree
    const unsigned int max_k = (tree.nparens / 2) - 1;

    const unsigned int num_prop_chunks = propstack_multi.get_num_stacks();
    while (k<max_k) {
          const unsigned int k_start = k;
          unsigned int filled_emb = 0;

          // chunk the progress to maximize cache reuse
#pragma omp parallel for 
          for (unsigned int ck=0; ck<num_prop_chunks; ck++) {
            su::PropStack<TFloat> &propstack = propstack_multi.get_prop_stack(ck);
            su::PropStack<TFloat> &countstack = countstack_multi.get_prop_stack(ck);
            const unsigned int tstart = propstack_multi.get_start(ck);
            const unsigned int tend = propstack_multi.get_end(ck);
            unsigned int my_filled_emb = 0;
            unsigned int my_k=k_start;

            while ((my_filled_emb<max_emb) && (my_k<max_k)) {
              const uint32_t node = tree.postorderselect(my_k);
              my_k++;

              TFloat *node_proportions = propstack.pop(node);
              TFloat *node_counts = countstack.pop(node);

              su::set_proportions_range(node_proportions, tree, node, table, tstart, tend, propstack);
              su::set_proportions_range(node_counts, tree, node, table, tstart, tend, countstack, false);

              if(task_p->bypass_tips && tree.isleaf(node))
                  continue;

              if (ck==0) { // they all do the same thing, so enough for the first to update the global state
                lengths[filled_emb] = tree.lengths[node];
                filled_emb++;
              }
              taskObj.embed_range(node_proportions, node_counts, tstart, tend, my_filled_emb);
              my_filled_emb++;
            }
            if (ck==0) { // they all do the same thing, so enough for the first to update the global state
              k=my_k;
            }
          }

#pragma acc wait
#pragma acc update device(lengths[:filled_emb])
          taskObj.sync_embedded(filled_emb);
          taskObj._run(filled_emb,lengths);
          filled_emb = 0;

          su::try_report(task_p, k, max_k);
    }

#pragma acc wait
    if(want_total) {
        const uint64_t start_idx = task_p->start;
        const uint64_t stop_idx = task_p->stop;

        TFloat * const dm_stripes_buf = taskObj.dm_stripes.buf;
        const TFloat * const dm_stripes_total_buf = taskObj.dm_stripes_total.buf;

#pragma acc parallel loop collapse(2) present(dm_stripes_buf,dm_stripes_total_buf)
        for(uint64_t i = start_idx; i < stop_idx; i++)
            for(uint64_t j = 0; j < n_samples; j++) {
                uint64_t idx = (i-start_idx)*n_samples_r+j;
                dm_stripes_buf[idx]=dm_stripes_buf[idx]/dm_stripes_total_buf[idx];
                // taskObj.dm_stripes[i][j] = taskObj.dm_stripes[i][j] / taskObj.dm_stripes_total[i][j];
            }

    }


#pragma acc exit data delete(lengths[:max_emb])
#pragma acc exit data delete(sample_total_counts[:n_samples_r])
    free(lengths);
    free(sample_total_counts);
}

void SUCMP_NM::unifrac_vaw(const su::biom_interface &table,
                           const su::BPTree &tree,
                           su::Method unifrac_method,
                           std::vector<double*> &dm_stripes,
                           std::vector<double*> &dm_stripes_total,
                           const su::task_parameters* task_p) {
    switch(unifrac_method) {
        case su::unweighted:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnweightedTask<double>,double>(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_normalized:
            unifrac_vawTT<SUCMP_NM::UnifracVawNormalizedWeightedTask<double>,double>(   table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_unnormalized:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnnormalizedWeightedTask<double>,double>( table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::generalized:
            unifrac_vawTT<SUCMP_NM::UnifracVawGeneralizedTask<double>,double>(          table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::unweighted_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnweightedTask<float >,float >(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_normalized_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawNormalizedWeightedTask<float >,float >(   table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_unnormalized_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnnormalizedWeightedTask<float >,float >( table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::generalized_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawGeneralizedTask<float >,float >(          table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        default:
            fprintf(stderr, "Unknown unifrac task\n");
            exit(1);
            break;
    }
}

