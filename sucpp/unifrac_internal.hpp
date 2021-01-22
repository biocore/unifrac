/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef __UNIFRAC_INTERNAL
#define __UNIFRAC_INTERNAL 1

#include <vector>
#include "biom_interface.hpp"
#include "task_parameters.hpp"

namespace su {
 // helper reporting function
 void try_report(const su::task_parameters* task_p, unsigned int k, unsigned int max_k);


 // Helper class with default constructor
 // The default is small enough to fit in L1 cache
 template<class TFloat>
 class PropStackFixed : public PropStack<TFloat> {
  public:
    static const uint32_t DEF_VEC_SIZE = 1024*sizeof(double)/sizeof(TFloat);

    PropStackFixed() : PropStack<TFloat>(DEF_VEC_SIZE) {}
 };

 // Helper class that splits a large vec_size into several smaller chunks of def_size
 template<class TFloat>
 class PropStackMulti {
  protected:
    const uint32_t vecsize;
    std::vector<PropStackFixed<TFloat> > multi;

  public:
    PropStackMulti(uint32_t _vecsize)
    : vecsize(_vecsize)
    , multi((vecsize + (PropStackFixed<TFloat>::DEF_VEC_SIZE-1))/PropStackFixed<TFloat>::DEF_VEC_SIZE) // round up
    {}
    ~PropStackMulti() {}

    uint32_t get_num_stacks() const {return (vecsize + (PropStackFixed<TFloat>::DEF_VEC_SIZE-1))/PropStackFixed<TFloat>::DEF_VEC_SIZE;}

    uint32_t get_start(uint32_t idx) const {return idx*PropStackFixed<TFloat>::DEF_VEC_SIZE;}
    uint32_t get_end(uint32_t idx) const   {return std::min((idx+1)*PropStackFixed<TFloat>::DEF_VEC_SIZE, vecsize);}

    PropStackFixed<TFloat> &get_prop_stack(uint32_t idx) {return multi[idx];}
 };

 void initialize_stripes(std::vector<double*> &dm_stripes,
                         std::vector<double*> &dm_stripes_total,
                         bool want_total,
                         const su::task_parameters* task_p);

 uint64_t initialize_sample_counts(double*& _counts, const su::task_parameters* task_p, const biom_interface &table);
 uint64_t initialize_sample_counts(float*& _counts, const su::task_parameters* task_p, const biom_interface &table);

}

#endif
