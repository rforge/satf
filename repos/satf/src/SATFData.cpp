#include <Rcpp.h>

using namespace Rcpp;

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include "SATFData.h"
#include "math_aux.h"


#define MAX_CORRELATION .96

#define GET(TYPE, NAME) Rcpp::clone( Rcpp::as<TYPE>(params[NAME]) )

#define debug_level 0

_dbg_class_set(CSATFResponseVariable, "CSATFResponseVariable", debug_level);
_dbg_class_set(CSATFDesignMatrix, "CSATFDesignMatrix", debug_level);
_dbg_class_set(CSATFData, "CSATFData", debug_level);



/**********************************
 *   class CSATFResponseVariable  *
 **********************************/

CSATFResponseVariable::CSATFResponseVariable(DoubleVector& predicted_criterion, CharacterVector& dv, 
                                             CharacterVector& cnames, DataFrame& data)
{
  _dbg_function("constructor", 1);

  mPredictedCriterion = Rcpp::clone( predicted_criterion );
  
  if( dv.containsElementNamed("response") ) {
    _dbg((1, "assigning 1"));
      std::string name_response = as< std::string >(dv["response"]);
      mResponseYes = as< std::vector<int> >(data[ name_response ]);
    _dbg((1, "/assigning 1"));
    mRVType = rv_binary;

  } else if( dv.containsElementNamed("dprime") ) {
      std::string name_dprime = as< std::string >(dv["dprime"]);
      mResponseDprime = as< std::vector<double> >(dv[ name_dprime ]);
      mRVType = rv_dprime;

  } else if( dv.containsElementNamed("n.responses") ) {
      std::string name_yes = as< std::string >(dv["n.responses.yes"]);
      std::string name_all = as< std::string >(dv["n.responses"]);
      mnResponsesYes = as< IntegerVector >(data[ name_yes ]);
      mnResponses = as< IntegerVector >(data[ name_all ]);
      mRVType = rv_aggregate;
  }
  
  if(cnames.containsElementNamed("trial.id") ) {
      std::string name_time = as< std::string >(cnames["trial.id"]);
      mTrialId = as< std::vector<int> >(data[name_time]);
  }
}

/**********************************
 *   class CDesignMatrix  *
 **********************************/


CSATFDesignMatrix::CSATFDesignMatrix(List& contrasts, DataFrame& data, 
                                     CharacterVector& cnames):  
                                        mCoefIndexFirst(parameter_invalid, 0),
                                        mCoefIndexLast(parameter_invalid, 0)
{
  _dbg_function("constructor", 1);

  std::string name_signal = as< std::string >(cnames["signal"]);
  IntegerVector signal = as< IntegerVector >( data[name_signal]);

  // determine the dimensions of the design matrix and set the indices
  CharacterVector names_contrast = contrasts.names();
  mnCoefs = 0;
  for(int i=0; i < contrasts.length(); i++) {
    std::string name = as< std::string >(names_contrast[i]);
    Parameter param_contrast = StringToParameter( name ); 
    CharacterVector terms = as< CharacterVector >(contrasts[i]);
    mCoefIndexFirst[param_contrast] = mnCoefs;
    mnCoefs += terms.length() + 1;
    mCoefIndexLast[param_contrast]  = mnCoefs-1;
  }
  mnDatapoints = data.nrows();

  // construct the design matrix
  NumericMatrix cur_dm = NumericMatrix(mnDatapoints, mnCoefs);
  int offset = 0;
  for(int i=0; i < contrasts.length(); i++) {
    std::string name_contrast = as< std::string >(names_contrast[i]);
    CharacterVector terms = contrasts[ name_contrast ];
    cur_dm(_, offset+0) = signal;
    for(int j=0; j < terms.length(); j++) {
      std::string cur_cname = as< std::string >(terms[j]);
      cur_dm(_, offset+j+1) = as< IntegerVector >( data[ cur_cname ] );
    }
    offset += 1 + terms.length();
  }
  mDM = cur_dm;
}

// void CSATFDesignMatrix::StringToParameter(std::string& name) {  
// }


CSATFDesignMatrix::Parameter CSATFDesignMatrix::StringToParameter(std::string& name)
{
  if(name == "lambda") {
    return parameter_lambda;

  } else if(name == "beta") {
    return parameter_beta;

  } else if(name == "delta") {
    return parameter_delta;

  }
  return parameter_invalid;
}

double CSATFDesignMatrix::GetParameter(Parameter parameter, int datapoint_index, std::vector<double>& coefs)  
{
    _dbg_function("GetParameter", 1);

    //_dbg((1, "nrow <%d>, ncol <%d>, coefslen: <%d>", mDM.nrow(), mDM.ncol(), coefs.size()));
    double param = 0.0;
    for(int idx = mCoefIndexFirst[parameter]; idx <= mCoefIndexLast[parameter]; idx++)
    {
      //_dbg((1, "row <%d>, col <%d>", datapoint_index, idx));
       double contrast = mDM(datapoint_index, idx);
      _dbg((0, "mDM(%d, %d) = %f", datapoint_index, idx, idx)); // TODO contrast));
       if(contrast != 0) {
          param += coefs[idx]*contrast;
       }
    }
    return(param);
}

double CSATFDesignMatrix::ComputeDprime(std::vector<double>& coefs, int datapoint_index, double time) 
{
  _dbg_function("ComputeDprime", 1);
  double lambda = GetParameter(parameter_lambda, datapoint_index, coefs);
  double beta   = GetParameter(parameter_beta,  datapoint_index, coefs);
  double delta  = GetParameter(parameter_delta, datapoint_index, coefs);
  _dbg((0, "time <%.2f>", time));
  _dbg((0, "coef <%.2f, %.2f, %.2f>", coefs[0], coefs[1], coefs[2]));
  _dbg((0, "lambda <%.2f>, beta <%.2f>, delta <%.2f>", lambda, beta, delta));
  return SATF(time, lambda, beta, delta);
}
  

/******************************
 *   class CSATFData          *
 ******************************/


CSATFData::CSATFData(CharacterVector& dv, List& contrasts,
                     NumericMatrix& coef_constraints, DataFrame& data, 
                     DoubleVector& predicted_criterion, 
                     CharacterVector& cnames):  
                            mDM(contrasts, data, cnames), 
                            mRV( predicted_criterion, dv, cnames, data ),
                            mDbgTime(20, 0)
{
  _dbg_function("constructor", 1);

    mParametersLower = coef_constraints(_, 0);
    mParametersUpper = coef_constraints(_, 1);
    
    std::string name_time = as< std::string >(cnames["time"]);
    mTime = as< std::vector<double> >(data[name_time]);
    
  //if( params.containsElementNamed("hyperparams") )
  //  mHyperparams = as<NumericMatrix>(params["hyperparams"]);
}

CSATFData::~CSATFData() {
//  _dbg_function("destructor", 1);
//  _dbg((-10, "mDbgTime[0]=%f sec", mDbgTime[0]));
}


inline double transform_with_one_boundary(double x, double lower, bool back) {
  if(!back) // [-Inf, +Inf] -> [lower, +Inf]
    return pow(x, 2)+lower;
  else // [lower, +Inf] -> [-Inf, +Inf]
    return sqrt(x - lower);
}

inline double transform_with_two_boundaries(double x, double lower, double upper, bool back) {
  if(!back) { // [-Inf, +Inf] -> [0, 1] -> [lower, upper]
    double y = exp(x)/(1+exp(x));
    double z = y*(upper-lower)+lower;
    return(z);
    
  } else {  // [lower, upper] -> [0, 1] -> [-Inf, +Inf]
    double y = (x-lower)/(upper-lower);
    double z = log(y/(1-y));
    return(z);
  }
}
double CSATFData::TransformCoef(double raw_coef, int i, bool back)
{
  _dbg_function("TransformCoef", 1);

  if(mParametersLower[i] == R_NegInf && mParametersUpper[i] == R_PosInf) {
    return raw_coef;     
  } 
  else if(mParametersUpper[i] == R_PosInf ) {
    return transform_with_one_boundary(raw_coef, mParametersLower[i], back);

  } else if( mParametersLower[i] == R_NegInf) {
    return -1*transform_with_one_boundary(raw_coef, -mParametersUpper[i], back);
    
  } else {
    return transform_with_two_boundaries(raw_coef, mParametersLower[i], mParametersUpper[i], back);
  }
}

std::vector<double> CSATFData::TransformCoefs(DoubleVector& raw_coefs, bool back)
{
  _dbg_function("TransformCoefs", 1);

  std::vector<double> coefs;
  int n_fixed = 0;
  for(int i=0; i < raw_coefs.size(); i++)
  {
    double val;
    if(mParametersLower[i] == mParametersUpper[i]) {
      val = mParametersLower[i];
      n_fixed++;
      
    } else {
      val = TransformCoef(raw_coefs[i-n_fixed], i, back);
    }
    coefs.push_back( val );  
  }
  return coefs;
}

double CSATFData::CoefsLL(std::vector<double>& coefs)
{
  _dbg_function("CoefsLL", 1);

  if( mHyperparams.length() == 0 )
    return 0;
    
  assert( mHyperparams.ncol() == 2 );
  assert( mHyperparams.nrow() == mDM.nrow() );

  double LL = 0;
  if( mHyperparams.nrow()*mHyperparams.ncol() == 0 )
    return 0;
  
  for(size_t i=0; i < mDM.nCoefs(); i++) {
    NumericVector hyperparams = mHyperparams( _, i);
    LL += _dgamma(coefs[i], hyperparams[0], hyperparams[1], true);
  }
  return LL;
}


double CSATFData::ObjectiveFunction_Binary(std::vector<double>& coefs)
{
  _dbg_function("ObjectiveFunction_Binary", 1);

  _dbg((0, "rv_type: %d, datapoints: %d", mRV.GetRVType(), mDM.nDatapoints()));
  assert( mRV.GetRVType() == rv_binary);

  double LL=0;
  double withinTrialCorrelation;
  double last_dprime = 0;
  
  if(coefs.size() > mDM.nCoefs()) {
    withinTrialCorrelation = logodds2p(coefs[coefs.size()]);
  }

  for(size_t idx = 0; idx < mDM.nDatapoints(); idx++)
  {
    _dbg((0, "time <%.1f>, crit <%.2f>, resp. <%d>", mTime[idx], mRV.PredictedCriterion(idx), mRV.ResponseYes(idx)));
    double dprime = mDM.ComputeDprime(coefs, idx, mTime[idx]);

    bool dependent = false;
    if(mRV.HasTrialId() && idx > 0 && withinTrialCorrelation != 0) {
      dependent = mRV.TrialId(idx)==mRV.TrialId(idx-1);
    }

    double cur_LL;
    if(dependent)
    {
      double CriterionMinusPsi = mRV.PredictedCriterion(idx) - dprime;
      double LastCriterionMinusPsi = mRV.PredictedCriterion(idx-1) - last_dprime;
  
      int lboundary_upper = 1 - mRV.ResponseYes(idx-1);
      double p_No = conditional_pnorm(CriterionMinusPsi, withinTrialCorrelation,
                                      LastCriterionMinusPsi, lboundary_upper);
      if(mRV.ResponseYes(idx) == 1) cur_LL = log(1-p_No);
      else                          cur_LL = log(p_No);
        
    } else 
    {
      int response = mRV.ResponseYes(idx);  
      double pNo = _pnorm(mRV.PredictedCriterion(idx), dprime);
      cur_LL = log( response - (2*response-1)*pNo );
     _dbg((0, "t <%.1f>, dprime <%.2f>, crit <%.2f>, resp. <%d>, LL <%.2f>, pNo <%.2f>", 
                  mTime[idx], dprime, mRV.PredictedCriterion(idx), mRV.ResponseYes(idx), cur_LL, pNo));
    }
    LL += cur_LL;
  }

  _dbg((0, "coefs <%.2f, %.2f, %.2f>", coefs[0], coefs[1], coefs[2])); 
  _dbg((0, "LL <%f>", LL)); 

  return LL;
}

double CSATFData::ObjectiveFunction_Aggregate(std::vector<double>& coefs)
{
  _dbg_function("ObjectiveFunction_Aggregate", 1);

  assert( mRV.GetRVType() == rv_aggregate);
  double LL = 0;
  
  for(size_t idx=0; idx < mDM.nDatapoints(); idx++)
  {
    double dprime = mDM.ComputeDprime(coefs, idx, mTime[idx]);
    
    double pYes = 1-_pnorm(mRV.PredictedCriterion(idx), dprime);
    LL += _dbinom(mRV.NResponsesYes(idx), mRV.NResponses(idx), pYes );
  }
  return LL;
}

double CSATFData::ObjectiveFunction(DoubleVector& raw_coefs)
{
  _dbg_function("ObjectiveFunction", 1);

  assert(coefs.length() == mnCoefs || coefs.length() == mnCoefs+1 );

  double LL = 0; 
  std::vector<double> coefs = TransformCoefs(raw_coefs);

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
