
#include "satf_math.h"

/* Approximation of the conditional bivariate normal cdf, based on Albers&Kallenberg (1994), 
   "A simple approximation to the bivariate normal distribution with large correlation coefficient".
    Other approximations are cited in Albers & Kallenberg (1994):
    * Cox & Wermuth (1991), and
    * http://www.dtic.mil/dtic/tr/fulltext/u2/a125033.pdf (when correlation is low)
    * http://www.nasa.gov/centers/dryden/pdf/87929main_H-1120.pdf 
    * Continuous Multivariate Distributions, Models and Applications (Kotz et al.)
*/

namespace AlbersKallenberg1994 {
  
  // derivative of the normal distribution density function
  inline double _ddnorm(double x) {
    return -_dnorm(x)*x;
  }
  inline double f_theta(double rho) {
    return( sqrt(1 - (rho*rho) )/rho );
  }

  inline double _pnorm2d(double a, double b, double rho, bool smaller, bool second_order)
  {
      double theta = sqrt(1 - pow(rho,2)) / rho;
      double c = (rho*a - b)/ sqrt(1-pow(rho,2));
      double c1 = _dnorm(c) / _pnorm(-c);
      double zeta = a + theta*(c1 - c);
      double eta = _pnorm(-c);

      double res;
      if(smaller) res = _pnorm(b);
      else        res = 1-_pnorm(a);
      res = res - eta*( _pnorm( zeta ) - _pnorm(a) );

      if(second_order)
        res = res - 0.5*(theta*theta)*eta*_ddnorm( zeta )*( 1 + (c - c1)*c1 );
      return res;
  }

  inline bool check_constraints(double a, double b, double rho) {
      assert( 0 < rho && rho < 1 );
      if( (rho*a-b) >= abs(rho*b - a)) return true;
      else                             return false;
  }
  
  inline double pnorm2d(double x_upper, double y_upper, double rho, bool second_order) {
      if( check_constraints(x_upper, y_upper, rho) ) {
        return _pnorm2d(x_upper, y_upper, rho, true, second_order);
      }
      
      if( check_constraints(y_upper, x_upper, rho) ){
        return _pnorm2d(y_upper, x_upper, rho, true, second_order);
      }
    
      if( check_constraints(-x_upper, -y_upper, rho) ) {
        return _pnorm2d(-x_upper, -y_upper, rho, false, second_order);
      }
    
      if( check_constraints(-y_upper, -x_upper, rho) ) {
        return _pnorm2d(-y_upper, -x_upper, rho, false, second_order);
      }
    
      return nan("");
  }
}

double pnorm2d(double x_upper, double y_upper, double rho, double second_order) {
    return AlbersKallenberg1994::pnorm2d(x_upper, y_upper, rho, second_order);
}

double pnorm_conditional(double rho, double relative_criterion, double last_relative_criterion, 
                         bool response_above_criterion, bool last_response_above_criterion, bool tolerate_imprecisions)
{
    if(last_response_above_criterion) {
      last_relative_criterion = -last_relative_criterion;
      relative_criterion = -relative_criterion;
    }
      
    double p_last = _pnorm(last_relative_criterion, 0.0, 1.0, true, false);
    double p_both = pnorm2d(relative_criterion, last_relative_criterion, rho, true);

    if(response_above_criterion == last_response_above_criterion) {
//      printf("%f(p_both)/%f(p_last) = %f\n", p_both, p_last, p_both/p_last);
      return log(p_both/p_last);
    }
//    printf("1 - %f(p_both)/%f(p_last) = %f\n", p_both, p_last, 1-p_both/p_last);
    return log(1 - p_both/p_last);
}

/*

P(X > x |Y > y) = P(X < -x |Y < -y)                     =     P(X < -x, Y < -y) / P(Y < -y)
P(X < x |Y > y) = P(X > -x |Y < -y) = P(X > -x |Y < -y) = 1 - P(X < -x, Y < -y) / P(Y < -y)
P(X > x |Y < y) = 1 - P(X < x |Y < y) = 1 - P(X < x, Y < y) / P(Y < y)
P(X < x |Y < y)                       =     P(X < x, Y < y) / P(Y < y)



*/

/*
  if( !valid_probability(pNo) ) {
    if(log_undefined_values) {
      printf("SATF WARNING: invalid conditional probability pNo=<%f>\n", pNo);
      printf("criterion=%f, dprime=%f, last_criterion=%f, last_dprime=%f\n", 
              criterion.value, dprime.value, last_datapoint.criterion.value, 
              last_datapoint.dprime.value);
      printf("corr.mrsat <%f>, relative_criterion <%.3f>, last_relative_criterion <%.3f>\n", corr_mrsat, relative_criterion, last_datapoint.relative_criterion);
      printf("last_response <%d>\n", last_datapoint.n_responses_yes);
    }
    return nan("");
  }
     

*/
