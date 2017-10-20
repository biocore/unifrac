#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "api.hpp"

void err(bool condition, const char* msg) {
    if(condition) {
        fprintf(stderr, "%s\n", msg);
        exit(1);
    }
} 


int main(int argc, char** argv) {
    mat_t* result = NULL;
    const char* table = "test.biom";
    const char* tree = "test.tre";
    const char* method = "unweighted";
    double exp[] = {0.2, 0.57142857, 0.6, 0.5, 0.2, 0.42857143, 0.66666667, 0.6, 0.33333333, 0.71428571, 0.85714286, 0.42857143, 0.33333333, 0.4, 0.6};
    
    ComputeStatus status;
    status = one_off(table, tree, method,
                     false, 1.0, false, 2, &result);

    err(status != okay, "Compute failed");
    err(result == NULL, "Empty result");
    err(result->n_samples != 6, "Wrong number of samples");
    err(result->cf_size != 15, "Wrong condensed form size");
    err(!result->is_upper_triangle, "Result is not squaure");

    for(unsigned int i = 0; i < result->cf_size; i++) 
        err(fabs(exp[i] - result->condensed_form[i]) > 0.00001, "Result is wrong");

    return 0;
}
