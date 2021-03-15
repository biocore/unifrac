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
#include "affinity.hpp"
#include <cstdlib>
#include <thread>
#include <signal.h>
#include <stdarg.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>

#include "unifrac_internal.hpp"

static pthread_mutex_t printf_mutex;
static bool* report_status;

static int sync_printf(const char *format, ...) {
    // https://stackoverflow.com/a/23587285/19741
    va_list args;
    va_start(args, format);

    pthread_mutex_lock(&printf_mutex);
    vprintf(format, args);
    pthread_mutex_unlock(&printf_mutex);

    va_end(args);
}

static void sig_handler(int signo) {
    // http://www.thegeekstuff.com/2012/03/catch-signals-sample-c-code
    if (signo == SIGUSR1) {
        if(report_status == NULL)
            fprintf(stderr, "Cannot report status.\n");
        else {
            for(int i = 0; i < CPU_SETSIZE; i++) {
                report_status[i] = true;
            }
        }
    }
}

using namespace su;

void su::try_report(const su::task_parameters* task_p, unsigned int k, unsigned int max_k) {
  if(__builtin_expect(report_status[task_p->tid], false)) {
    sync_printf("tid:%u\tstart:%u\tstop:%u\tk:%u\ttotal:%u\n", task_p->tid, task_p->start, task_p->stop, k, max_k);
    report_status[task_p->tid] = false;
  }
}

void su::register_report_status() {
    // register a signal handler so we can ask the master thread for its
    // progress
    if (signal(SIGUSR1, sig_handler) == SIG_ERR)
        fprintf(stderr, "Can't catch SIGUSR1\n");

    report_status = (bool*)calloc(sizeof(bool), CPU_SETSIZE);
    pthread_mutex_init(&printf_mutex, NULL);
}

void su::remove_report_status() {
    if(report_status != NULL) {
        pthread_mutex_destroy(&printf_mutex);
        free(report_status);
        report_status = NULL;
    }
}

template<class TFloat>
PropStack<TFloat>::PropStack(uint32_t vecsize) 
: prop_stack()
, prop_map()
, defaultsize(vecsize)
{
    prop_map.reserve(1000);
}

template<class TFloat>
PropStack<TFloat>::~PropStack() {
    // drain stack
    for(unsigned int i = 0; i < prop_stack.size(); i++) {
        TFloat *vec = prop_stack.top();
        prop_stack.pop();
        free(vec);
    }

    // drain the map
    for(auto it = prop_map.begin(); it != prop_map.end(); it++) {
        TFloat *vec = it->second;
        free(vec);
    }
    prop_map.clear();
}

template<class TFloat>
TFloat* PropStack<TFloat>::get(uint32_t i) {
    return prop_map[i];
}

template<class TFloat>
void PropStack<TFloat>::push(uint32_t node) {
    TFloat* vec = prop_map[node];
    prop_map.erase(node);
    prop_stack.push(vec);
}

template<class TFloat>
TFloat* PropStack<TFloat>::pop(uint32_t node) {
    /*
     * if we don't have any available vectors, create one
     * add it to our record of known vectors so we can track our mallocs
     */
    TFloat *vec;
    int err = 0;
    if(prop_stack.empty()) {
        err = posix_memalign((void **)&vec, 32, sizeof(TFloat) * defaultsize);
        if(vec == NULL || err != 0) {
            fprintf(stderr, "Failed to allocate %zd bytes, err %d; [%s]:%d\n",
                    sizeof(TFloat) * defaultsize, err, __FILE__, __LINE__);
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

// make sure they get instantiated
template class su::PropStack<float>;
template class su::PropStack<double>;


void su::initialize_stripes(std::vector<double*> &dm_stripes,
                            std::vector<double*> &dm_stripes_total,
                            bool want_total,
                            const su::task_parameters* task_p) {
    int err = 0;
    for(unsigned int i = task_p->start; i < task_p->stop; i++){
        err = posix_memalign((void **)&dm_stripes[i], 4096, sizeof(double) * task_p->n_samples);
        if(dm_stripes[i] == NULL || err != 0) {
            fprintf(stderr, "Failed to allocate %zd bytes, err %d; [%s]:%d\n",
                    sizeof(double) * task_p->n_samples, err, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        for(unsigned int j = 0; j < task_p->n_samples; j++)
            dm_stripes[i][j] = 0.;

        if(want_total) {
            err = posix_memalign((void **)&dm_stripes_total[i], 4096, sizeof(double) * task_p->n_samples);
            if(dm_stripes_total[i] == NULL || err != 0) {
                fprintf(stderr, "Failed to allocate %zd bytes err %d; [%s]:%d\n",
                        sizeof(double) * task_p->n_samples, err, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            for(unsigned int j = 0; j < task_p->n_samples; j++)
                dm_stripes_total[i][j] = 0.;
        }
    }
}

template<class TFloat>
void su::set_proportions(TFloat* __restrict__ props,
                         const BPTree &tree,
                         uint32_t node,
                         const biom_interface &table,
                         PropStack<TFloat> &ps,
                         bool normalize) {
    if(tree.isleaf(node)) {
       table.get_obs_data(tree.names[node], props);
       if (normalize) {
#pragma omp parallel for schedule(static)
        for(unsigned int i = 0; i < table.n_samples; i++) {
           props[i] /= table.sample_counts[i];
        }
       }

    } else {
        unsigned int current = tree.leftchild(node);
        unsigned int right = tree.rightchild(node);

#pragma omp parallel for schedule(static)
        for(unsigned int i = 0; i < table.n_samples; i++)
            props[i] = 0;

        while(current <= right && current != 0) {
            TFloat * __restrict__ vec = ps.get(current);  // pull from prop map
            ps.push(current);  // remove from prop map, place back on stack

#pragma omp parallel for schedule(static)
            for(unsigned int i = 0; i < table.n_samples; i++)
                props[i] = props[i] + vec[i];

            current = tree.rightsibling(current);
        }
    }
}

// make sure they get instantiated
template void su::set_proportions(float* __restrict__ props,
                                  const BPTree &tree,
                                  uint32_t node,
                                  const biom_interface &table,
                                  PropStack<float> &ps,
                                  bool normalize);
template void su::set_proportions(double* __restrict__ props,
                                  const BPTree &tree,
                                  uint32_t node,
                                  const biom_interface &table,
                                  PropStack<double> &ps,
                                  bool normalize);

template<class TFloat>
void su::set_proportions_range(TFloat* __restrict__ props,
                               const BPTree &tree,
                               uint32_t node,
                               const biom_interface &table, 
                               unsigned int start, unsigned int end,
                               PropStack<TFloat> &ps,
                               bool normalize) {
    const unsigned int els = end-start;
    if(tree.isleaf(node)) {
       table.get_obs_data_range(tree.names[node], start, end, normalize, props);
    } else {
        const unsigned int right = tree.rightchild(node);
        unsigned int current = tree.leftchild(node);

        for(unsigned int i = 0; i < els; i++)
            props[i] = 0;

        while(current <= right && current != 0) {
            const TFloat * __restrict__ vec = ps.get(current);  // pull from prop map
            ps.push(current);  // remove from prop map, place back on stack

            for(unsigned int i = 0; i < els; i++)
                props[i] += vec[i];

            current = tree.rightsibling(current);
        }
    }
}

// make sure they get instantiated
template void su::set_proportions_range(float* __restrict__ props,
                                        const BPTree &tree,
                                        uint32_t node,
                                        const biom_interface &table,
                                        unsigned int start, unsigned int end,
                                        PropStack<float> &ps,
                                        bool normalize);
template void su::set_proportions_range(double* __restrict__ props,
                                        const BPTree &tree,
                                        uint32_t node,
                                        const biom_interface &table,
                                        unsigned int start, unsigned int end,
                                        PropStack<double> &ps,
                                        bool normalize);

std::vector<double*> su::make_strides(unsigned int n_samples) {
    uint32_t n_rotations = (n_samples + 1) / 2;
    std::vector<double*> dm_stripes(n_rotations);

    int err = 0;
    for(unsigned int i = 0; i < n_rotations; i++) {
        double* tmp;
        err = posix_memalign((void **)&tmp, 32, sizeof(double) * n_samples);
        if(tmp == NULL || err != 0) {
            fprintf(stderr, "Failed to allocate %zd bytes, err %d; [%s]:%d\n",
                    sizeof(double) * n_samples, err,  __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        for(unsigned int j = 0; j < n_samples; j++)
            tmp[j] = 0.0;
        dm_stripes[i] = tmp;
    }
    return dm_stripes;
}

