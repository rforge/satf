
#include <Rcpp.h>
using namespace Rcpp;
#include <stdio.h>
#include "satf.h"
#include "satf_math.h"

#include "debug.h"

int global_dbg_level = 10;

static CSATFData *zzz=NULL;

// [[Rcpp::export]]
LogicalVector rcpp_initialize_logLikFn(CharacterVector& dv, NumericMatrix& dm, IntegerVector& dm_ncoef, 
                                       NumericMatrix& constraints, DataFrame& data, CharacterVector& cnames)
{	
  _dbg_function("rcpp_initialize_logLikFn", 1);
    zzz = new CSATFData(dv, dm, dm_ncoef, constraints, data, cnames);
  _dbg((0, "returning"));
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
LogicalVector rcpp_deinitialize_logLikFn() {
  _dbg_function("rcpp_deinitialize_logLikFn", 1);
  _dbg((0, "pointer <0x%x>", zzz));
  	delete zzz;
  	zzz = NULL;
  _dbg((0, "returning"));
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_logLikFn(DoubleVector& coefs, bool by_row=false, bool tolerate_imprecision=true) {
  _dbg_function("rcpp_compute_logLikFn", 1);
  
  // DoubleVector gradients(coefs.size());
  DoubleVector res = zzz->ObjectiveFunction(coefs, by_row, tolerate_imprecision, NULL);
  // res.attr("gradient") = gradients;
  _dbg((0, "returning"));
  return Rcpp::wrap( res );
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_logLikFn_gradient(DoubleVector& coefs, bool by_row=false, bool tolerate_imprecision=true) {
  _dbg_function("rcpp_compute_logLikFn_gradient", 1);
  DoubleVector gradients(coefs.size());
  DoubleVector LL = zzz->ObjectiveFunction(coefs, by_row, tolerate_imprecision, &gradients);
  _dbg((0, "returning"));
  return Rcpp::wrap(gradients);
}

// [[Rcpp::export]]
DoubleVector rcpp_constrain_coefs(DoubleVector& coefs) {
  _dbg_function("rcpp_constrain_coefs", 1);
  _dbg((0, "init"));
  DoubleVector res = zzz->ConstrainCoefs(coefs, true);
  _dbg((0, "returning"));
  return Rcpp::wrap(res);
}

// [[Rcpp::export]]
DoubleVector rcpp_unconstrain_coefs(DoubleVector& coefs) {
  _dbg_function("rcpp_constrain_coefs", 1);
  _dbg((0, "init"));
  DoubleVector res = zzz->UnconstrainCoefs(coefs);
  _dbg((0, "returning"));
  return Rcpp::wrap(res);
}

// [[Rcpp::export]]
void rcpp_select_subset_by_zero_dm_columns_any(LogicalVector& zero_column) {
  _dbg_function("rcpp_select_subset_by_zero_dm_columns_any", 1);
  zzz->SelectSubset(zero_column, false);
  _dbg((0, "returning"));
}

// [[Rcpp::export]]
void rcpp_select_subset_by_zero_dm_columns_all(LogicalVector& zero_column) {
  _dbg_function("rcpp_select_subset_by_zero_dm_columns_any", 1);
  zzz->SelectSubset(zero_column, true);
  _dbg((0, "returning"));
}

// [[Rcpp::export]]
void rcpp_reset_selection( ) {
  _dbg_function("rcpp_reset_selection", 1);
  zzz->ResetSubset();
  _dbg((0, "returning"));
}

// [[Rcpp::export]]
void rcpp_set_coef_values(DoubleVector& values) {
  _dbg_function("rcpp_set_coef_values", 1);
  zzz->SetCoefValues(values);
  _dbg((0, "returning"));
}

// [[Rcpp::export]]
void rcpp_reset_coef_values( CharacterVector& names ) {
  _dbg_function("rcpp_reset_coef_values", 1);
   zzz->ResetCoefValues(names);
  _dbg((0, "returning"));
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

/*
// [[Rcpp::export]]
NumericMatrix rcpp_get_dm() {
  return Rcpp::wrap(zzz->mDM.mDM);
}

// [[Rcpp::export]]
NumericVector rcpp_get_constraints_lower() {
  return Rcpp::wrap(zzz->mCoefConstraints.mCoefsLower);
}
// [[Rcpp::export]]
NumericVector rcpp_get_constraints_upper() {
  return Rcpp::wrap(zzz->mCoefConstraints.mCoefsUpper);
}
*/

#if 0
  #include "satf.cpp"
  #include "satf_math.cpp"
  #include "debug.cpp"
#endif
