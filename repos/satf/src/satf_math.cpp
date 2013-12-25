
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
  inline double _ddnorm(double x, double mu, double sigma) {
    return -_dnorm(x, mu, sigma)*(x-mu)/pow(sigma,2);
  }

  inline double _ddnorm(double x) {
    return -_dnorm(x)*x;
  }

  inline double f_c(double a, double b, double rho) {
    return (rho*a - b)/sqrt(1-(rho*rho));
  }
  inline double f_c1(double c) {
    return _dnorm(c) / _pnorm(-c);
  }
  inline double f_theta(double rho) {
    return( sqrt(1-pow(rho,2))/rho );
  }
  inline double g(double a, double b, double rho, double c, double c1, double theta) {
    return( _pnorm(-c)*( _pnorm(a + theta*(c1 - c)) - _pnorm(a) ) ); 
  }
  inline double h(double a, double b, double rho, double c, double c1, double theta) {
    return 0.5*pow(theta,2)*_pnorm(-c)*_ddnorm(a+theta*(c1-c))*(1+c*c1-pow(c1,2));
  }

  inline bool check_constraints(double a, double b, double rho) {
      assert( 0 < rho && rho < 1 );
      if( (rho*a-b) >= abs(rho*b - a)) return true;
      else                             return false;
  }

  inline double _pnorm2d(double a, double b, double rho, double theta, 
                         bool smaller, bool second_order)
  {
      if( !check_constraints(a, b, rho) )
          return(-1);
          
      double res;
      if(smaller) res = _pnorm(b);
      else        res = (1-_pnorm(a));
          
      double c = AlbersKallenberg1994::f_c(a,b,rho);
      double c1 = AlbersKallenberg1994::f_c1(c);
      res = res - g(a, b, rho, c, c1, theta);
      if(second_order)
        return res - h(a, b, rho, c, c1, theta);
      else
        return res;
  }
  
  inline double pnorm2d(double x_upper, double y_upper, double rho, bool second_order) {
      double theta, res;
      theta = AlbersKallenberg1994::f_theta(rho);
      
      res = _pnorm2d(x_upper, y_upper, rho, theta, true, second_order);
      if(res >= 0) return res;
      
      res = _pnorm2d(y_upper, x_upper, rho, theta, true, second_order);
      if(res >= 0) return res;
    
      res = _pnorm2d(-x_upper, -y_upper, rho, theta, false, second_order);
      if(res >= 0) return res;
    
      res = _pnorm2d(-y_upper, -x_upper, rho, theta, false, second_order);
      if(res >= 0) return res;
    
      return nan("");
  }
}

double pnorm2d(double x_upper, double y_upper, double rho, double second_order) {
    return AlbersKallenberg1994::pnorm2d(x_upper, y_upper, rho, second_order);
}

double pnorm_conditional(double rho, double crit_minus_psi, double last_crit_minus_psi, 
                         bool last_response_above_criterion, bool tolerate_imprecisions)
{
    double p_below;
    if(last_response_above_criterion) {
        p_below = 1-pnorm2d(-crit_minus_psi, -last_crit_minus_psi, rho)/_pnorm(-last_crit_minus_psi);
        
    } else {
        p_below = pnorm2d(crit_minus_psi, last_crit_minus_psi, rho)/_pnorm(last_crit_minus_psi);
    }

    if(tolerate_imprecisions) {
        static const double max_error = 0.0005;
        static const double default_prob = 0.0001;
        if(p_below < 0 && p_below > -max_error) {
          p_below = default_prob;
          
        } else if(p_below > 1 && p_below < 1+max_error) {
          p_below = 1-default_prob;
        }
    }    
    return p_below;
}
