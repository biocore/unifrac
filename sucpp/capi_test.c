#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "api.hpp"

#ifndef bool
#define bool char
#define true 1
#define false 0
#endif

void err(bool condition, const char* msg) {
    if(condition) {
        fprintf(stderr, "%s\n", msg);
        exit(1);
    }
} 

void test_su(int num_cores){
    mat_t* result = NULL;
    const char* table = "test.biom";
    const char* tree = "test.tre";
    const char* method = "unweighted";
    double exp[] = {0.2, 0.57142857, 0.6, 0.5, 0.2, 0.42857143, 0.66666667, 0.6, 0.33333333, 0.71428571, 0.85714286, 0.42857143, 0.33333333, 0.4, 0.6};
    
    ComputeStatus status;
    status = one_off(table, tree, method,
                     false, 1.0, false, num_cores, &result);

    err(status != okay, "Compute failed");
    err(result == NULL, "Empty result");
    err(result->n_samples != 6, "Wrong number of samples");
    err(result->cf_size != 15, "Wrong condensed form size");
    err(!result->is_upper_triangle, "Result is not squaure");

    for(unsigned int i = 0; i < result->cf_size; i++) 
        err(fabs(exp[i] - result->condensed_form[i]) > 0.00001, "Result is wrong");

}

void test_faith_pd(){
    r_vec* result = NULL;
    const char* table = "test.biom";
    const char* tree = "test.tre";
    double exp[] = {4, 5, 6, 3, 2, 5};

    ComputeStatus status;
    status = faith_pd_one_off(table, tree, &result);

    err(status != okay, "Compute failed");
    err(result == NULL, "Empty result");
    err(result->n_samples != 6, "Wrong number of samples");

    for(unsigned int i = 0; i < result->n_samples; i++)
        err(fabs(exp[i] - result->values[i]) > 0.00001, "Result is wrong");

}

int main(int argc, char** argv) {
    int num_cores = strtol(argv[1], NULL, 10);

    printf("Testing Striped UniFrac...\n");
    test_su(num_cores);
    printf("Tests passed.\n");
    printf("Testing Faith's PD...\n");
    test_faith_pd();
    printf("Tests passed.\n");
    return 0;
}

