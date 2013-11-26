
#include <Rcpp.h>
using namespace Rcpp;
#include <stdio.h>
#include "SATFData.h"
#include "satf_math.h"

static CSATFData *zzz=NULL;

// [[Rcpp::export]]
LogicalVector rcpp_initialize_logLikFn(CharacterVector& dv,
					List& contrasts, 
					NumericMatrix& coef_constraints, DataFrame& data, 
					DoubleVector& predicted_criterion, 
					CharacterVector& cnames) {	
    zzz = new CSATFData(dv, contrasts, coef_constraints, data, 
					              predicted_criterion, cnames);
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
LogicalVector rcpp_deinitialize_logLikFn() {
	delete zzz;
	zzz = NULL;
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_logLikFn(DoubleVector& coefs, bool by_row=false) {
    return Rcpp::wrap(zzz->ObjectiveFunction(coefs, by_row) );
}

// [[Rcpp::export]]
DoubleVector rcpp_constrain_coefs(DoubleVector& coefs) {
    return Rcpp::wrap(zzz->ConstrainCoefs(coefs, true));
}

// [[Rcpp::export]]
DoubleVector rcpp_unconstrain_coefs(DoubleVector& coefs) {
    return Rcpp::wrap(zzz->UnconstrainCoefs(coefs));
}

// [[Rcpp::export]]
DoubleVector rcpp_pnorm2d(DoubleVector& x_lower, DoubleVector& y_lower, DoubleVector& rho, LogicalVector& order) {
    return Rcpp::wrap(pnorm2d(x_lower[0], y_lower[0], rho[0], order[0]));
}

// [[Rcpp::export]]
DoubleVector rcpp_correlate(IntegerVector& trial_id, DoubleVector& noise, DoubleVector& rho_vec) {
    double rho = rho_vec[0];
    for(int i=1; i < trial_id.length(); i++) {
      if(trial_id[i-1] == trial_id[i])
        noise[i] = noise[i-1]*rho + noise[i]*sqrt(1-pow(rho,2));
    }
  return Rcpp::wrap(noise);
}

/*
// [[Rcpp::export]]
void rcpp_add_coefficient(DoubleVector& lower, DoubleVector& upper, CharacterVector& name) {
  std::string name_str = as<std::string>(name[0]);
  zzz->AddCoefficient((double)lower[0], (double)upper[0], name_str );
}
*/

#if 1
#include "satf.cpp"
#include "satf_math.cpp"
#
