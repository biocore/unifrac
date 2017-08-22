#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include <fstream>
#include <string>
#include <iomanip>
#include <thread>
#include <vector>

#define OK 0
#define TREE_MISSING 1
#define TABLE_MISSING 2
#define UNKNOWN_METHOD 3

namespace su {
    /* a result matrix
     *
     * n_samples <uint> the number of samples.
     * cf_size <uint> the size of the condensed form.
     * is_square <bool> if true, indicates condensed_form represents a square
     *      matrix, and only the upper triangle is contained. if false, 
     *      condensed_form represents a dissimilarity matrix in which case 
     *      condensed_form stores all non-diagonal values in row order.
     * condensed_form <double*> the matrix values of length cf_size.
     * sample_ids <char**> the sample IDs of length n_samples.
     */
    /*struct mat {
        unsigned int n_samples;  
        unsigned int cf_size;   
        bool is_square;          
        double* condensed_form;
        char** sample_ids;
    };*/
    

    int one_off(const char* biom_filename, const char* tree_filename, 
                const char* unifrac_method, bool variance_adjust, double alpha,
                unsigned int threads, double** &result);
}

