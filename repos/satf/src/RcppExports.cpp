// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_initialize_logLikFn
LogicalVector rcpp_initialize_logLikFn(CharacterVector& dv, NumericMatrix& dm, IntegerVector& dm_ncoef, NumericMatrix& constraints, DataFrame& data, CharacterVector& cnames);
RcppExport SEXP satf_rcpp_initialize_logLikFn(SEXP dvSEXP, SEXP dmSEXP, SEXP dm_ncoefSEXP, SEXP constraintsSEXP, SEXP dataSEXP, SEXP cnamesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< CharacterVector& >::type dv(dvSEXP );
        Rcpp::traits::input_parameter< NumericMatrix& >::type dm(dmSEXP );
        Rcpp::traits::input_parameter< IntegerVector& >::type dm_ncoef(dm_ncoefSEXP );
        Rcpp::traits::input_parameter< NumericMatrix& >::type constraints(constraintsSEXP );
        Rcpp::traits::input_parameter< DataFrame& >::type data(dataSEXP );
        Rcpp::traits::input_parameter< CharacterVector& >::type cnames(cnamesSEXP );
        LogicalVector __result = rcpp_initialize_logLikFn(dv, dm, dm_ncoef, constraints, data, cnames);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_deinitialize_logLikFn
LogicalVector rcpp_deinitialize_logLikFn();
RcppExport SEXP satf_rcpp_deinitialize_logLikFn() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        LogicalVector __result = rcpp_deinitialize_logLikFn();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_compute_logLikFn
DoubleVector rcpp_compute_logLikFn(DoubleVector& coefs, bool by_row = false, bool tolerate_imprecision = true);
RcppExport SEXP satf_rcpp_compute_logLikFn(SEXP coefsSEXP, SEXP by_rowSEXP, SEXP tolerate_imprecisionSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< DoubleVector& >::type coefs(coefsSEXP );
        Rcpp::traits::input_parameter< bool >::type by_row(by_rowSEXP );
        Rcpp::traits::input_parameter< bool >::type tolerate_imprecision(tolerate_imprecisionSEXP );
        DoubleVector __result = rcpp_compute_logLikFn(coefs, by_row, tolerate_imprecision);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_compute_logLikFn_gradient
DoubleVector rcpp_compute_logLikFn_gradient(DoubleVector& coefs, bool by_row = false, bool tolerate_imprecision = true);
RcppExport SEXP satf_rcpp_compute_logLikFn_gradient(SEXP coefsSEXP, SEXP by_rowSEXP, SEXP tolerate_imprecisionSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< DoubleVector& >::type coefs(coefsSEXP );
        Rcpp::traits::input_parameter< bool >::type by_row(by_rowSEXP );
        Rcpp::traits::input_parameter< bool >::type tolerate_imprecision(tolerate_imprecisionSEXP );
        DoubleVector __result = rcpp_compute_logLikFn_gradient(coefs, by_row, tolerate_imprecision);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_constrain_coefs
DoubleVector rcpp_constrain_coefs(DoubleVector& coefs);
RcppExport SEXP satf_rcpp_constrain_coefs(SEXP coefsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< DoubleVector& >::type coefs(coefsSEXP );
        DoubleVector __result = rcpp_constrain_coefs(coefs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_unconstrain_coefs
DoubleVector rcpp_unconstrain_coefs(DoubleVector& coefs);
RcppExport SEXP satf_rcpp_unconstrain_coefs(SEXP coefsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< DoubleVector& >::type coefs(coefsSEXP );
        DoubleVector __result = rcpp_unconstrain_coefs(coefs);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_select_subset_by_zero_dm_columns_any
void rcpp_select_subset_by_zero_dm_columns_any(LogicalVector& zero_column);
RcppExport SEXP satf_rcpp_select_subset_by_zero_dm_columns_any(SEXP zero_columnSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< LogicalVector& >::type zero_column(zero_columnSEXP );
        rcpp_select_subset_by_zero_dm_columns_any(zero_column);
    }
    return R_NilValue;
END_RCPP
}
// rcpp_select_subset_by_zero_dm_columns_all
void rcpp_select_subset_by_zero_dm_columns_all(LogicalVector& zero_column);
RcppExport SEXP satf_rcpp_select_subset_by_zero_dm_columns_all(SEXP zero_columnSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< LogicalVector& >::type zero_column(zero_columnSEXP );
        rcpp_select_subset_by_zero_dm_columns_all(zero_column);
    }
    return R_NilValue;
END_RCPP
}
// rcpp_reset_selection
void rcpp_reset_selection();
RcppExport SEXP satf_rcpp_reset_selection() {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        rcpp_reset_selection();
    }
    return R_NilValue;
END_RCPP
}
// rcpp_set_coef_values
void rcpp_set_coef_values(DoubleVector& values);
RcppExport SEXP satf_rcpp_set_coef_values(SEXP valuesSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< DoubleVector& >::type values(valuesSEXP );
        rcpp_set_coef_values(values);
    }
    return R_NilValue;
END_RCPP
}
// rcpp_reset_coef_values
void rcpp_reset_coef_values(CharacterVector& names);
RcppExport SEXP satf_rcpp_reset_coef_values(SEXP namesSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< CharacterVector& >::type names(namesSEXP );
        rcpp_reset_coef_values(names);
    }
    return R_NilValue;
END_RCPP
}
// rcpp_pnorm2d
double rcpp_pnorm2d(double x_lower, double y_lower, double rho, bool second_order);
RcppExport SEXP satf_rcpp_pnorm2d(SEXP x_lowerSEXP, SEXP y_lowerSEXP, SEXP rhoSEXP, SEXP second_orderSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type x_lower(x_lowerSEXP );
        Rcpp::traits::input_parameter< double >::type y_lower(y_lowerSEXP );
        Rcpp::traits::input_parameter< double >::type rho(rhoSEXP );
        Rcpp::traits::input_parameter< bool >::type second_order(second_orderSEXP );
        double __result = rcpp_pnorm2d(x_lower, y_lower, rho, second_order);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_correlate
DoubleVector rcpp_correlate(IntegerVector& trial_id, DoubleVector& noise, DoubleVector& rho_vec);
RcppExport SEXP satf_rcpp_correlate(SEXP trial_idSEXP, SEXP noiseSEXP, SEXP rho_vecSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector& >::type trial_id(trial_idSEXP );
        Rcpp::traits::input_parameter< DoubleVector& >::type noise(noiseSEXP );
        Rcpp::traits::input_parameter< DoubleVector& >::type rho_vec(rho_vecSEXP );
        DoubleVector __result = rcpp_correlate(trial_id, noise, rho_vec);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
