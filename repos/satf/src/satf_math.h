#ifndef __MATH_AUX_H__
#define __MATH_AUX_H__

#include <Rcpp.h>
#include <assert.h>

double pnorm_conditional(double rho, double relative_criterion, double last_relative_criterion, 
                         bool response_above_criterion, bool last_response_above_criterion, bool tolerate_imprecisions=true);
double pnorm2d(double x_upper, double y_upper, double rho, double second_order=true);

inline double logodds2p(double lodds) { return( exp(lodds)/(1+exp(lodds)) ); }
inline double _dnorm(double x, double mu=0.0, double sigma=1.0, bool lg=false) { return ::Rf_dnorm4(x, mu, sigma, lg?1:0); }
inline double _pnorm(double x, double mu=0.0, double sigma=1.0, bool lt=true, int lg=false) { return ::Rf_pnorm5(x, mu, sigma, lt?1:0, lg?1:0); }
inline double _dgamma(double x, double shp, double scl, bool lg=false) { return ::Rf_dgamma(x, shp, scl, lg?1:0); }
inline double _dbinom(double x, double n, double p, bool lg=false)     { return ::Rf_dbinom(x, n, p, lg?1:0); }


// negatively accelerated exponential
inline double NAE(double t, double asymptote, double invrate, double intercept) {
  if(asymptote == 0)
    return 0.0;
  else if(t <= intercept)
    return 0.0;
  else if(invrate <= 0)
    return nan("");
  else
    return asymptote*(1-exp(-(1/invrate)*(t-intercept)));
}

// shifted negatively accelerated exponential
inline double SNAE(double time, double max, double invrate, double intercept, double min) {
	return NAE(time, max-min, invrate, intercept) + min;
}


class SNAEPoint {
public:
  inline SNAEPoint() {
    reset();
  }
  
  inline SNAEPoint(double t, double fasymptote, double finvrate, double fintercept, double fmin=0.0) {
    Update(t, fasymptote, finvrate, fintercept, fmin);
  }
  
  inline void Update(double t, double fasymptote, double finvrate, double fintercept, double fmin=0.0) {
    time = t;
    asymptote = fasymptote;
    invrate = finvrate;
    intercept = fintercept;
    min = fmin;
    if(asymptote == fmin)
      value = fmin;
    else
      value = SNAE(time, asymptote, invrate, intercept, min);
  }
  inline void reset() {
    time = nan("");
    asymptote = nan("");
    min = nan("");
    invrate = nan("");
    intercept = nan("");
    value = nan("");
  }
  
public:
  double time;
  double asymptote;
  double min;
  double invrate;
  double intercept;
  double value;
};


#endif // __MATH_AUX_H__
