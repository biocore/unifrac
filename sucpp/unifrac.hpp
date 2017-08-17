#include <stack>
#include <vector>
#include <unordered_map>
#include <thread>
#include "unifrac_task.hpp"

namespace su {
    enum Method {unweighted, weighted_normalized, weighted_unnormalized, generalized};
    
    class PropStack {
        private:
            std::stack<double*> prop_stack;
            std::unordered_map<uint32_t, double*> prop_map;
            uint32_t defaultsize;
        public:
            PropStack(uint32_t vecsize);
            ~PropStack();
            double* pop(uint32_t i);
            void push(uint32_t i);
            double* get(uint32_t i);
    };

    void unifrac(biom &table, 
                 BPTree &tree, 
                 Method unifrac_method,
                 std::vector<double*> &dm_stripes,
                 std::vector<double*> &dm_stripes_total,
                 const task_parameters* task_p);
    
    void unifrac_vaw(biom &table, 
                     BPTree &tree, 
                     Method unifrac_method,
                     std::vector<double*> &dm_stripes,
                     std::vector<double*> &dm_stripes_total,
                     const task_parameters* task_p);
    
    double** deconvolute_stripes(std::vector<double*> &stripes, uint32_t n);
    void set_proportions(double* props, 
                         BPTree &tree, uint32_t node, 
                         biom &table, 
                         PropStack &ps,
                         bool normalize = true);
    std::vector<double*> make_strides(unsigned int n_samples);
    inline void embed_proportions(double* out, double* in, uint32_t n) {
        double val;
        for(unsigned int i = 0; i < n; i++) {
            val = in[i];
            out[i] = val;
            out[i + n] = val;
        }
    }

}
