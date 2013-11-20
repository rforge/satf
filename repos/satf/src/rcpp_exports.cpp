
#include <Rcpp.h>
using namespace Rcpp;
#include <stdio.h>
#include "SATFData.h"

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
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
DoubleVector rcpp_compute_logLikFn(DoubleVector& coefs) {
/*	clock_t t = clock();

	DoubleVector x;
	for(int i=0; i < 200; i++) {
		x = zzz->ObjectiveFunction(coefs);
	}

	double ft = ((float)(clock()-t))/CLOCKS_PER_SEC;
	printf("ticks <%f> sec\n\n", ft);
*/
    return Rcpp::wrap(zzz->ObjectiveFunction(coefs) );
}

// [[Rcpp::export]]
DoubleVector rcpp_transform_coefs(DoubleVector& coefs) {
    return Rcpp::wrap(zzz->TransformCoefs(coefs, false));
}

// [[Rcpp::export]]
DoubleVector rcpp_untransform_coefs(DoubleVector& coefs) {
    return Rcpp::wrap(zzz->TransformCoefs(coefs, true));
}

// [[Rcpp::export]]
List rcpp_xxx() {
    return Rcpp::wrap( List::create( zzz->mDM.mDM ) );
}

#include "SATFData.cpp"