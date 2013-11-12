
#include <Rcpp.h>
using namespace Rcpp;
#include <stdio.h>
#include "SATFData.h"

static CSATFData *zzz;

// [[Rcpp::export]]
LogicalVector rcpp_initialize_logLikFn(List& params) {
    zzz = new CSATFData(params);
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
LogicalVector rcpp_deinitialize_logLikFn() {
    delete zzz;
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_logLikFn(DoubleVector& coefs, DoubleVector& fixed_coefs) {
    return Rcpp::wrap(zzz->ObjectiveFunction(coefs, fixed_coefs) );
}

//#include "SATFData.cpp"