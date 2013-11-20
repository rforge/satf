#ifndef __MATH_AUX_H__
#define __MATH_AUX_H__

#include <Rcpp.h>

inline double logodds2p(double lodds) { return( exp(lodds)/(1+exp(lodds)) ); }
inline double _dnorm(double x, double mu=0.0, double sigma=1.0, bool lg=false) { return ::Rf_dnorm4(x, mu, sigma, lg?1:0); }
inline double _pnorm(double x, double mu=0.0, double sigma=1.0, bool lt=true, int lg=false) { return ::Rf_pnorm5(x, mu, sigma, lt?1:0, lg?1:0); }
inline double _dgamma(double x, double shp, double scl, bool lg=false) { return ::Rf_dgamma(x, shp, scl, lg?1:0); }
inline double _dbinom(double x, double n, double p, bool lg=false)     { return ::Rf_dbinom(x, n, p, lg?1:0); }

// inline int index(int row, int col, int nRows) {
//    return( (col*nRows + row)-1 );
// }

inline double SATF(double t, double lambda, double beta, double delta) {
  if(lambda == 0)
    return 0.0;
  return( (t >= delta)*(lambda*(1-exp(-(1/beta)*(t-delta)))) );
}


/* Approximation of the conditional bivariate normal density based on the following
   paper: http://www.dtic.mil/dtic/tr/fulltext/u2/a125033.pdf (approximation to X, 
   given Y < b, best when correlation is low)
   Background: http://www.nasa.gov/centers/dryden/pdf/87929main_H-1120.pdf (density of X, 
   given Y < b).  Auch: Continuous Multivariate Distributions, Models and Applications (Kotz et al.)
*/


inline double conditional_pnorm(double nboundary, double rho,
                                double lboundary, int lboundary_upper)
{
  int mirror = (lboundary_upper*2-1);
  rho = rho*mirror;
  lboundary = lboundary*mirror;
  double mu = -rho*_dnorm(lboundary, 0, 1, false)/_pnorm(lboundary);
  double sigma = sqrt(1 + rho*lboundary*mu - pow(mu, 2) );
  return( _pnorm((nboundary-mu)/sigma) );
}


#endif // __MATH_AUX_H__
