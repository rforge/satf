#include <Rcpp.h>

using namespace Rcpp;

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <vector>

#include "SATFData.h"
#include "math_aux.h"


#define MAX_CORRELATION .96

#define DEBUG printf

#define GET(TYPE, NAME) Rcpp::clone( Rcpp::as<TYPE>(params[NAME]) )
    
CSATFResponseVariable::CSATFResponseVariable(List& allparams)
{
  mPredictedCriterion = Rcpp::clone( Rcpp::as<DoubleVector>(allparams["predicted.criterion"]) );    
  
  List params = Rcpp::as<List>(allparams["rv"]);

  mRVType = (enumRVType) Rcpp::as<int>(params["rv.type"]);
  switch(mRVType) 
  {
    case rv_binary:
      mResponseYes = GET(IntegerVector, "response.yes");
      mTrialId = GET(IntegerVector, "trial.id");
      break;
      
    case rv_aggregate:
      mResponseDprime = GET(DoubleVector, "response.dprime");
      break;
      
    case rv_dprime:
      mResponseHits = GET(DoubleVector, "n.response.hits");
      mnGrammatical = GET(DoubleVector, "n.grammatical");
      mResponseFAs = GET(DoubleVector, "n.response.fas");
      mnUngrammatical = GET(DoubleVector, "n.ungrammatical");
      break;
  }
}


// TODO: check that all other parameters agree with these numbers of coefficients and data points
// TODO: add checks in the wrapper functions which make sure the asserts here don't fail

CSATFData::CSATFData( List& params ):  mDM( params ), mRV( params )
{
    mFixedparams = GET(DoubleVector, "fixedparams");

    mTransformVec = GET(IntegerVector, "transform.vec");
    mParamsMinimum = GET(DoubleVector, "params.min");
    mParamsMaximum = GET(DoubleVector, "params.max");
        
    mTime = GET(DoubleVector, "time");
    
    mHyperparams = GET(NumericMatrix, "hyperparams");
    
    assert(mTransformVec.length() == mnCoefs);
    assert(mParamsMinimum.length() == mnCoefs);
    assert(mParamsMaximum.length() == mnCoefs);
}



CSATFDesignMatrix::CSATFDesignMatrix( Rcpp::List& allparams )
{
    List params = Rcpp::as<List>(allparams["dm"]);

    mDM = GET(NumericMatrix, "design.matrix");
    
    mRowsLambda = GET(IntegerVector, "lambda");
    mRowsBeta = GET(IntegerVector, "beta");
    mRowsDelta = GET(IntegerVector, "delta");
    
    assert(mTransformVec.length() == mnCoefs);
    assert(mParamsMinimum.length() == mnCoefs);
    assert(mParamsMaximum.length() == mnCoefs);
}


double CSATFDesignMatrix::GetParam(Rcpp::IntegerVector& rows, int datapoint_index, Rcpp::DoubleVector& coefs)  
{
   double param = 0.0;
   for(int coef_idx=0; coef_idx < rows.length(); coef_idx++) {
       int row = rows[coef_idx]-1;
       double contrast = mDM(row, datapoint_index);
       if(contrast) {
          param += coefs[row]*contrast;
       }
   }
   return(param);
}

double CSATFDesignMatrix::ComputeDprime(DoubleVector& coefs, int datapoint_index, double time) {
    double lambda = GetParam(mRowsLambda, datapoint_index, coefs);
    double beta = GetParam(mRowsBeta,  datapoint_index, coefs);
    double delta = GetParam(mRowsDelta, datapoint_index, coefs);
    return SATF(time, lambda, beta, delta);
}
  



// TODO: Fix the interaction of maximim and minimum values with transformations here: Currently, it is a hack.

void CSATFData::CoefsTransform(DoubleVector& coefs, IntegerVector& transform_vec, 
  				                     DoubleVector& minimum, DoubleVector& maximum)
{
  int n_proper_coefs = transform_vec.length();

  for(int i=0; i < n_proper_coefs; i++) {
    if(transform_vec[i] != 0) {
      coefs[i] = pow(coefs[i],2)*transform_vec[i];
    }
    if(coefs[i] < minimum[i]) {
      coefs[i] = minimum[i];
    }
     else if(coefs[i] > maximum[i]) {
      coefs[i] = maximum[i];
    }
  }

  if( coefs.length() > n_proper_coefs ) {
    coefs[n_proper_coefs] = logodds2p(coefs[n_proper_coefs]);
    if(coefs[n_proper_coefs] > MAX_CORRELATION) {
      coefs[n_proper_coefs] = MAX_CORRELATION;
    }
  }
}




double CSATFData::CoefsLL(DoubleVector& coefs)
{
  assert( mHyperparams.ncol() == 2 );
  assert( mHyperparams.nrow() == mDM.nrow() );

  double LL = 0;
  if( mHyperparams.nrow()*mHyperparams.ncol() == 0 )
    return 0;
  
  for(int i=0; i < mnCoefs; i++) {
    NumericVector hyperparams = mHyperparams( _, i);
    LL += _dgamma(coefs[i], hyperparams[0], hyperparams[1], true);
  }
  return LL;
}


double CSATFData::ObjectiveFunction_Binary(DoubleVector& coefs)
{
  assert( mRV.GetRVType() == rv_binary);

  double LL=0;
  double withinTrialCorrelation;
  double last_dprime = 0;

  if(coefs.length() > mnCoefs) {
    withinTrialCorrelation = logodds2p(coefs[coefs.length()]);
  }
  
  for(int idx=0; idx < mnDatapoints; idx++)
  {
    double dprime = mDM.ComputeDprime(coefs, idx, mTime[idx]);
  
    bool dependent = FALSE;
    if(mRV.HasTrialId() && idx > 0 && withinTrialCorrelation != 0) {
      dependent = mRV.TrialId(idx)==mRV.TrialId(idx-1);
    }
    
    if(dependent)
    {
      double CriterionMinusPsi = mRV.PredictedCriterion(idx) - dprime;
      double LastCriterionMinusPsi = mRV.PredictedCriterion(idx-1) - last_dprime;
  
      int lboundary_upper = 1 - mRV.ResponseYes(idx-1);
      double p_No = conditional_pnorm(CriterionMinusPsi, withinTrialCorrelation,
                                      LastCriterionMinusPsi, lboundary_upper);
      if(mRV.ResponseYes(idx) == 1) LL += log(1-p_No);
      else                          LL += log(p_No);
        
    } else {
      int response = mRV.ResponseYes(idx);  
      double pNo = _pnorm(mRV.PredictedCriterion(idx), dprime);
      LL += log( response - (2*response-1)*pNo );
    }
  }
  return LL;
}

double CSATFData::ObjectiveFunction_Aggregate(DoubleVector& coefs)
{
  assert( mRV.GetRVType() == rv_aggregate);
  double LL = 0;
  
  for(int idx=0; idx < mnDatapoints; idx++)
  {
    double dprime = mDM.ComputeDprime(coefs, idx, mTime[idx]);
    
    double pYes = 1-_pnorm(mRV.PredictedCriterion(idx), dprime);
    LL += _dbinom(mRV.ResponseHits(idx), mRV.NGrammatical(idx), pYes );
    if(mRV.NUngrammatical(idx) != 0 ) 
    {
      double pYes = 1-_pnorm(mRV.PredictedCriterion(idx));
      LL += _dbinom(mRV.ResponseFAs(idx), mRV.NGrammatical(idx), pYes );
    }
  }
  return LL;
}

double CSATFData::ObjectiveFunction(DoubleVector& coefs, Rcpp::DoubleVector& fixed_coefs)
{
  assert(coefs.length() == mnCoefs || coefs.length() == mnCoefs+1 );

  double LL = 0; 

  CoefsTransform(coefs, mTransformVec, mParamsMinimum, mParamsMaximum);

  double LLparams = CoefsLL(coefs);
  LL += LLparams;
  
  switch(mRV.GetRVType()) 
  {
      case rv_aggregate:
        LL += ObjectiveFunction_Aggregate(coefs);
        break;
        
      case rv_binary:
        LL += ObjectiveFunction_Binary(coefs);
        break;
        
      case rv_dprime:
        // TODO: Not implemented for now.
        break;
  }  
  return(LL);
}
