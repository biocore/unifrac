#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "api.hpp"

void err(const char* msg) {
	fprintf(stderr, "%s\n", msg);
	exit(1);
} 


int main(int argc, char** argv) {
	mat_t* result = NULL;
	const char* table = "test.biom";
	const char* tree = "test.tre";
	const char* method = "unweighted";
	double exp[] = {0.2, 0.57142857, 0.6, 0.5, 0.2, 0.42857143, 0.66666667, 0.6, 0.33333333, 0.71428571, 0.85714286, 0.42857143, 0.33333333, 0.4, 0.6};
    
	ComputeStatus status;
    status = one_off(table, tree, method,
                     false, 1.0, 2, &result);

	if(status != okay)
		err("Compute failed");
	
	if(result == NULL) 
		err("Empty result");

	if(result->n_samples != 6) 
		err("Wrong number of samples");

	if(result->cf_size != 15) 
		err("Wrong condensed form size");
   
	if(!result->is_square)
		err("Result is not squaure");

	for(unsigned int i = 0; i < result->cf_size; i++) {
		if(fabs(exp[i] - result->condensed_form[i]) > 0.00001) {
			err("Result is wrong");
		}
	}

	return 0;
}
