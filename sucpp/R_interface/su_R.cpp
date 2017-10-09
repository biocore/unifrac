#include <iostream>
#include <Rcpp.h>
#include <vector>
#include "../api.hpp"

using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List unifrac(const char* table, const char* tree, int nthreads){
    mat_t* result = NULL;
    const char* method = "unweighted";
    ComputeStatus status;
    status = one_off(table, tree, method, false, 1.0, nthreads, &result);
    vector<double> cf;
    //push result->condensed_form into a vector becuase R doesn't like double*
    for(int i=0; i<result->cf_size; i++){
        cf.push_back(result->condensed_form[i]);
    }

    return Rcpp::List::create(Rcpp::Named("n_samples") = result->n_samples,
                              Rcpp::Named("is_sqaure") = result->is_square,
                              Rcpp::Named("cf_size") = result->cf_size,
                              Rcpp::Named("c_form") = cf);

}


