
#include <Rcpp.h>
#include "satf.h"
using namespace Rcpp;


CDesignMatrixRow::CDesignMatrixRow(DoubleVector& elements, 
                                   CCoefConstraints& coef_constraints):
                                   mCurrentParameters(parameter_invalid, nan(""))
{
  _dbg_method("constructor", 1);

  mElements = elements;
  mCoefConstraints = &coef_constraints;
}

bool CDesignMatrixRow::UpdateParameters(CCoefs& coefs, bool force_update)
{
  _dbg_method("UpdateParameters", 1);
  bool params_changed = false;
  for(int p = 0; p < parameter_invalid; p++) {
    if( force_update || coefs.HasParameterChanged((Parameter)p, *this) ) {
      mCurrentParameters[p] = coefs.ComputeParameter((Parameter)p, *this);
      params_changed = true;
    }
  }
  return params_changed;
}

double CDesignMatrixRow::operator==(CDesignMatrixRow& other) const
{ 
  if(mElements.size() != other.mElements.size())
    return false;
    
  for(int i=0; i < other.mElements.size(); i++) {
    if( mElements[i] != other.mElements[i] )
      return false;
  }
  return true;
}

