#ifndef __MATH_AUX_H__
#define __MATH_AUX_H__

#include <Rcpp.h>
#include <assert.h>

double pnorm_conditional(double rho, double crit_minus_psi, double last_crit_minus_psi, bool last_response_above_criterion);
double pnorm2d(double x_upper, double y_upper, double rho, double second_order=true);

inline double logodds2p(double lodds) { return( exp(lodds)/(1+exp(lodds)) ); }
inline double _dnorm(double x, double mu=0.0, double sigma=1.0, bool lg=false) { return ::Rf_dnorm4(x, mu, sigma, lg?1:0); }
inline double _pnorm(double x, double mu=0.0, double sigma=1.0, bool lt=true, int lg=false) { return ::Rf_pnorm5(x, mu, sigma, lt?1:0, lg?1:0); }
inline double _dgamma(double x, double shp, double scl, bool lg=false) { return ::Rf_dgamma(x, shp, scl, lg?1:0); }
inline double _dbinom(double x, double n, double p, bool lg=false)     { return ::Rf_dbinom(x, n, p, lg?1:0); }

inline double SATF(double t, double asymptote, double invrate, double intercept) {
  if(asymptote == 0)
    return 0.0;
  else if(t <= intercept)
    return 0.0;
  else 
    return asymptote*(1-exp(-(1/invrate)*(t-intercept)));
}

inline double CriterionF(double t, double upper, double lower, double center, double stretch) {
  return exp( (t-center)/stretch )/(1+exp( (t-center)/stretch ))*(upper-lower)+lower;
}

#endif // __MATH_AUX_H__
