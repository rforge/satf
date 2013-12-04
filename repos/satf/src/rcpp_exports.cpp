
#include <Rcpp.h>
using namespace Rcpp;
#include <stdio.h>
#include "satf.h"
#include "satf_math.h"

static CSATFData *zzz=NULL;

// [[Rcpp::export]]
LogicalVector rcpp_initialize_logLikFn(CharacterVector& dv, NumericMatrix& dm,
          IntegerVector& dm_ncoef, NumericMatrix& constraints, DataFrame& data, 
					CharacterVector& cnames)
{	
    zzz = new CSATFData(dv, dm, dm_ncoef, constraints, data, cnames);
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
LogicalVector rcpp_deinitialize_logLikFn() {
  	delete zzz;
  	zzz = NULL;
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_logLikFn(DoubleVector& coefs, bool by_row=false, bool tolerate_imprecision=true) {
  return Rcpp::wrap(zzz->ObjectiveFunction(coefs, by_row, tolerate_imprecision) );
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
void rcpp_select_subset(LogicalVector& selection) {
  zzz->SelectSubset(selection);
}

// [[Rcpp::export]]
void rcpp_reset_selection( ) {
  zzz->ResetSubset();
}

// [[Rcpp::export]]
void rcpp_set_coef_values(DoubleVector& values) {
  zzz->SetCoefValues(values);
}

// [[Rcpp::export]]
void rcpp_reset_coef_values( CharacterVector& names ) {
    zzz->ResetCoefValues(names);
}


// [[Rcpp::export]]
double rcpp_pnorm2d(double x_lower, double y_lower, double rho, bool second_order) {
    return pnorm2d(x_lower, y_lower, rho, second_order);
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

// [[Rcpp::export]]
NumericMatrix rcpp_get_dm() {
  return Rcpp::wrap(zzz->mDM.mDM);
}


/*
// [[Rcpp::export]]
void rcpp_add_coefficient(DoubleVector& lower, DoubleVector& upper, CharacterVector& name) {
  std::string name_str = as<std::string>(name[0]);
  zzz->AddCoefficient((double)lower[0], (double)upper[0], name_str );
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_criterion_logLikFn(DoubleVector& coefs, DoubleVector& time, 
                                             IntegerVector& response, IntegerVector& trial_id) {
  return Rcpp::wrap(ObjectiveFunctionCriterion(coefs, time, response, trial_id));
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_criterion(DoubleVector& coefs, DoubleVector& time) {
  return Rcpp::wrap(ComputeCriterion(coefs, time));
}
*/



#if 0
  #include "satf.cpp"
  #include "satf_math.cpp"
  #include "debug.cpp"
#endif
