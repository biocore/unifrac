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
#include <stack>
#include <unordered_map>
#include "biom_interface.hpp"
#include "task_parameters.hpp"
#include "unifrac.hpp"

namespace su {
 // helper reporting functions
 void register_report_status();
 void remove_report_status();
 void try_report(const su::task_parameters* task_p, unsigned int k, unsigned int max_k);

 template<class TFloat>
 class PropStack {
   private:
     std::stack<TFloat*> prop_stack;
     std::unordered_map<uint32_t, TFloat*> prop_map;
     uint32_t defaultsize;
   public:
     PropStack(uint32_t vecsize);
     virtual ~PropStack();
     TFloat* pop(uint32_t i);
     void push(uint32_t i);
     TFloat* get(uint32_t i);
 };

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

 template<class TFloat>
 void set_proportions(TFloat* __restrict__ props,
                      const BPTree &tree, uint32_t node,
                      const biom_interface &table,
                      PropStack<TFloat> &ps,
                      bool normalize = true);

 template<class TFloat>
 void set_proportions_range(TFloat* __restrict__ props,
                            const BPTree &tree, uint32_t node,
                            const biom_interface &table,unsigned int start, unsigned int end,
                            PropStack<TFloat> &ps,
                            bool normalize = true);


 void initialize_stripes(std::vector<double*> &dm_stripes,
                         std::vector<double*> &dm_stripes_total,
                         bool want_total,
                         const su::task_parameters* task_p);

  std::vector<double*> make_strides(unsigned int n_samples);

}

#endif
