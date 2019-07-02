#include <iostream>
#include <Rcpp.h>
#include <vector>
#include "api.hpp"

using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List unifrac(const char* table, const char* tree, int nthreads){
    mat_t* result = NULL;
    const char* method = "unweighted";
    ComputeStatus status;
    status = one_off(table, tree, method, false, 1.0, false, nthreads, &result);
    vector<double> cf;
    //push result->condensed_form into a vector becuase R doesn't like double*
    for(int i=0; i<result->cf_size; i++){
        cf.push_back(result->condensed_form[i]);
    }

    return Rcpp::List::create(Rcpp::Named("n_samples") = result->n_samples,
                              Rcpp::Named("is_upper_triangle") = result->is_upper_triangle,
                              Rcpp::Named("cf_size") = result->cf_size,
                              Rcpp::Named("c_form") = cf);

}

// [[Rcpp::export]]
Rcpp::List faith_pd(const char* table, const char* tree){
    r_vec* result = NULL;
    ComputeStatus status;
    status = faith_pd_one_off(table, tree, &result);
    vector<double> values;
    for(int i = 0; i < result->n_samples; i++){
        values.push_back(result->values[i]);
    }

    return Rcpp::List::create(Rcpp::Named("n_samples") = result->n_samples,
                              Rcpp::Named("faith_pd") = values);

}


