#include <Rcpp.h>

using namespace Rcpp;

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <sys/timeb.h>
#include <vector>

#include "satf.h"
#include "satf_math.h"



#define debug_level 10

_dbg_class_set(CCoefConstraints, "CCoefConstraints", 10);
_dbg_class_set(CDependentVariable, "CDependentVariable", 0);
_dbg_class_set(CDesignMatrix, "CDesignMatrix", 0);
_dbg_class_set(CSATFData, "CSATFData", 10);


int getMilliCount(){
  timeb tb;
	ftime(&tb);
	int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
	return nCount;
}

int getMilliSpan(int nTimeStart){
	int nSpan = getMilliCount() - nTimeStart;
	if(nSpan < 0)
		nSpan += 0x100000 * 1000;
	return nSpan;
}


/**********************************
 *   class CDependentVariable  *
 **********************************/

CDependentVariable::CDependentVariable(CharacterVector& dv, CharacterVector& cnames, DataFrame& data)
{
  _dbg_method("constructor", 1);

  if( dv.containsElementNamed("response") ) {
    _dbg((1, "assigning 1"));
      std::string name_response = as< std::string >(dv["response"]);
      mResponseYes = as< std::vector<int> >(data[ name_response ]);
    _dbg((1, "/assigning 1"));
    mDVType = rv_binary;

  } else if( dv.containsElementNamed("dprime") ) {
      std::string name_dprime = as< std::string >(dv["dprime"]);
      mResponseDprime = as< std::vector<double> >(dv[ name_dprime ]);
      mDVType = rv_dprime;

  } else if( dv.containsElementNamed("n.responses") ) {
      std::string name_yes = as< std::string >(dv["n.responses.yes"]);
      std::string name_all = as< std::string >(dv["n.responses"]);
      mnResponsesYes = as< IntegerVector >(data[ name_yes ]);
      mnResponses = as< IntegerVector >(data[ name_all ]);
      mDVType = rv_aggregate;
  }
  
  if(cnames.containsElementNamed("trial.id") ) {
      std::string name_time = as< std::string >(cnames["trial.id"]);
      mTrialId = as< std::vector<int> >(data[name_time]);
  }
}

CDependentVariable::~CDependentVariable()
{
  _dbg_method("destructor", 1);
}

/**********************************
 *   class CDesignMatrix  *
 **********************************/

CDesignMatrix::CDesignMatrix(Rcpp::NumericMatrix& dm, Rcpp::IntegerVector& dm_ncoef):  
                                        mCoefTypes(dm.ncol(), -1),
                                        mCoefIndexFirst(parameter_invalid, -1),
                                        mCoefIndexLast(parameter_invalid, -1)
{
  _dbg_method("constructor", 1);

  // copy design matrix
  mDM = Rcpp::clone(dm);


  // create a mapping from parameter type to coefficient indices
  CharacterVector ncoef_names = dm_ncoef.names();
  int cur_coef_index = 0;
  for(int i=0; i < dm_ncoef.size(); i++) {
    std::string cur_name = as<std::string>(ncoef_names[i]);
    Parameter parameter = StringToParameter( cur_name );
    mCoefIndexFirst[parameter] = cur_coef_index;
    mCoefIndexLast[parameter] = mCoefIndexFirst[parameter] + dm_ncoef[i] - 1;
    cur_coef_index += dm_ncoef[i];
  }

  // create a mapping from coefficient index to parameter type
  for(size_t i=0; i < mCoefIndexFirst.size(); i++) {
      if(mCoefIndexFirst[i] < 0)
        continue;
        
      for(int j=mCoefIndexFirst[i]; j <= mCoefIndexLast[i]; j++) {
          mCoefTypes[j] = i;
      }
  }

  // save coefficient names
  List dimnames(dm.attr("dimnames"));
  if(dimnames.size() > 1) 
  {
    RObject names = dimnames[1];
    if (!names.isNULL()) {
      CharacterVector coef_names = CharacterVector(names);
      _dbg((0, "copying <%d> names", coef_names.size() ));
      for(int i=0; i < coef_names.size(); i++) {
        mCoefNames.push_back( as< std::string >(coef_names[i]) );
      }
    }
  }

}


CDesignMatrix::~CDesignMatrix()
{
  _dbg_method("destructor", 1);
  _dbg((0, "destructing"));
}


CDesignMatrix::Parameter CDesignMatrix::StringToParameter(std::string& name)
{
  if(name == "asymptote") {
    return parameter_satf_asymptote;

  } else if(name == "invrate") {
    return parameter_satf_invrate;

  } else if(name == "intercept") {
    return parameter_satf_intercept;

  } else if(name == "bias.min") {
    return parameter_bias_min;

  } else if(name == "bias.max") {
    return parameter_bias_max;

  } else if(name == "bias.invrate") {
    return parameter_bias_invrate;
    
  } else if(name == "bias.intercept") {
    return parameter_bias_intercept;
    
  } else if(name == "corr.mrsat") {
    return parameter_corr_mrsat;
  }
  return parameter_invalid;
}

double CDesignMatrix::GetParameter(Parameter parameter, CCoefs& coefs, int datapoint_index)  
{
    _dbg_method("GetParameter", 1);
    double cur_coef, param = 0.0;

    if(!HasParameter(parameter))
      return nan("");

    //_dbg((1, "nrow <%d>, ncol <%d>, coefslen: <%d>", mDM.nrow(), mDM.ncol(), coefs.size()));
    for(int idx = mCoefIndexFirst[parameter]; idx <= mCoefIndexLast[parameter]; idx++)
    {
      //_dbg((1, "row <%d>, col <%d>", datapoint_index, idx));
       cur_coef = coefs[idx];
       if(cur_coef == 0.0)
          continue;
       double contrast = mDM(datapoint_index, idx);
       _dbg((0, "mDM(%d, %d) = %f", datapoint_index, idx, contrast));
       param += cur_coef*contrast;
    }
    return(param);
}

double CDesignMatrix::ComputeDprime(CCoefs& coefs, int datapoint_index, double time, bool log) 
{
  _dbg_method("ComputeDprime", 1);
  double asymptote = GetParameter(parameter_satf_asymptote, coefs, datapoint_index);
  if(asymptote == 0.0)
    return 0.0;
  if(isnan(asymptote))
    return asymptote;
  double invrate   = GetParameter(parameter_satf_invrate,  coefs, datapoint_index);
  double intercept = GetParameter(parameter_satf_intercept, coefs, datapoint_index);
  if(log)
    printf("SATF(asymptote=%f, invrate=%f, intercept=%f) = %f\n", asymptote, invrate, intercept, NAE(time, asymptote, invrate, intercept) );
  return NAE(time, asymptote, invrate, intercept);
}

// This is what Wickens (2001) calls lambda_center in equation (2.5) 
double CDesignMatrix::ComputeCriterion(CCoefs& coefs, int datapoint_index, double time, bool log) 
{
  _dbg_method("ComputeCprime", 1);
  double min = GetParameter(parameter_bias_min, coefs, datapoint_index);
  double max = GetParameter(parameter_bias_max, coefs, datapoint_index);
  double invrate = GetParameter(parameter_bias_invrate, coefs, datapoint_index);
  double intercept = GetParameter(parameter_bias_intercept, coefs, datapoint_index);
  if(log)
    printf("criterion(max=%f, invrate=%f, intercept=%f, min=%f) = %f\n", max, invrate, intercept, min, SNAE(time, max, invrate, intercept, min) );
  return SNAE(time, max, invrate, intercept, min);
}


double CDesignMatrix::ComputeLogLik_Independent(CCoefs& coefs, int datapoint_index, double time, 
                                                int n_responses_yes, int n_responses, DoubleVector* gradient) 
{
	_dbg_method("ComputeLogLik_Independent", 1);

  double dprime_asymptote = GetParameter(parameter_satf_asymptote, coefs, datapoint_index);
  double dprime_invrate   = GetParameter(parameter_satf_invrate,  coefs, datapoint_index);
  double dprime_intercept = GetParameter(parameter_satf_intercept, coefs, datapoint_index);
  double dprime = NAE(time, dprime_asymptote, dprime_invrate, dprime_intercept);

  double criterion_max = GetParameter(parameter_bias_max, coefs, datapoint_index);
  double criterion_invrate = GetParameter(parameter_bias_invrate, coefs, datapoint_index);
  double criterion_intercept = GetParameter(parameter_bias_intercept, coefs, datapoint_index);
  double criterion_min = GetParameter(parameter_bias_min, coefs, datapoint_index);
  double criterion = SNAE(time, criterion_max, criterion_invrate, criterion_intercept, criterion_min);

  if(gradient != NULL) {
    ComputeGradient(coefs, datapoint_index, dprime_asymptote, dprime_invrate, dprime_intercept, 
                    criterion_max, criterion_invrate, criterion_intercept, criterion_min,
                    criterion, dprime, time, n_responses_yes, n_responses, *gradient);
  }

	double pYes = _pnorm(-(criterion-dprime) );
	return _dbinom(n_responses_yes, n_responses, pYes, true );      
}


void CDesignMatrix::ComputeGradient(CCoefs& coefs, int datapoint_index,
                                     double dprime_asymptote, double dprime_invrate, double dprime_intercept, 
                                     double criterion_max, double criterion_invrate, double criterion_intercept, double criterion_min,
                                     double criterion, double dprime, double t, int n_responses_yes, int n_responses, DoubleVector& gradients)
{
    if( gradients.size() != coefs.size() )
      return;

    double criterion_tau = exp(-(t-criterion_intercept)/criterion_invrate);
    double criterion_at = (criterion_max- criterion_min)*criterion_tau - criterion_max + dprime; 
    double criterion_dnorm_at = _dnorm(criterion_at);
    double criterion_pnorm_at = _pnorm(criterion_at);
    double criterion_eta = criterion_pnorm_at*(-1+criterion_pnorm_at);

    double dprime_tau = exp(-(t-dprime_intercept)/dprime_invrate);
    double dprime_at = -criterion + dprime_asymptote - dprime_asymptote*dprime_tau;
    double dprime_dnorm_at = _dnorm(dprime_at);
    double dprime_pnorm_at = _pnorm(dprime_at);
    double dprime_eta = dprime_pnorm_at * (-1 + dprime_pnorm_at);

    for(int coef_idx = 0; coef_idx < coefs.size(); coef_idx++)
    {
        double gradient = 0;

        double w = mDM(datapoint_index, coef_idx);
        double criterion_zeta = w * criterion_dnorm_at * criterion_tau;
        double dprime_zeta = w * dprime_dnorm_at * dprime_asymptote * dprime_tau;

        if(w != 0)
        {
          switch(mCoefTypes[coef_idx])
          {
              case parameter_satf_asymptote:
                if(t < dprime_intercept)
                  continue;

                gradient = -(dprime_dnorm_at * w) / dprime_eta;
                gradient = gradient * (n_responses_yes*(1- dprime_tau) 
                                    + (-1 + dprime_tau)*n_responses * dprime_pnorm_at );
                break;

              case parameter_satf_invrate:
                if(t < dprime_intercept)
                  continue;

                gradient = -dprime_zeta / ( dprime_eta * pow(dprime_invrate,2) );
                gradient = gradient * ((-t+dprime_intercept) * n_responses_yes 
                                       + n_responses * dprime_pnorm_at * (t-dprime_intercept));
                break;

              case parameter_satf_intercept:
                if(t < dprime_intercept)
                  continue;

                gradient = -dprime_zeta / ( dprime_eta * dprime_invrate );
                gradient = gradient * (-n_responses_yes + n_responses * dprime_pnorm_at);
                break;

              case parameter_bias_max:
                if(t < criterion_intercept)
                  continue;

                gradient = criterion_dnorm_at * w / criterion_eta;
                gradient = gradient * (n_responses_yes - n_responses*criterion_pnorm_at  + 
                                      (n_responses*criterion_pnorm_at - n_responses_yes)*criterion_tau );
                break;
 
            case parameter_bias_invrate:
              if(t < criterion_intercept)
                continue;
              {
                double t_minus_intercept = (t - criterion_intercept);
                double intercept_minus_t = (criterion_intercept - t);

                gradient = criterion_zeta / (criterion_eta * pow(criterion_invrate,2));
                gradient = gradient*(n_responses_yes*(criterion_max*intercept_minus_t + criterion_min*t_minus_intercept)
                            + n_responses*criterion_pnorm_at*(criterion_max*t_minus_intercept + criterion_min*intercept_minus_t));
              }
                break;

            case parameter_bias_intercept:
              if(t < criterion_intercept)
                continue;
               gradient = criterion_zeta / ( criterion_eta * criterion_invrate );
               gradient = gradient*(n_responses_yes*(-criterion_max+criterion_min) + n_responses*criterion_pnorm_at*(criterion_max-criterion_min));
              break;

            case parameter_bias_min:
              if(t > criterion_intercept) {
                  gradient = -criterion_zeta / criterion_eta;
                  gradient = gradient*( -n_responses_yes + n_responses*criterion_pnorm_at );

              } else  {
                  double dnorm_at = _dnorm(-criterion + dprime);
                  double pnorm_at = _pnorm(-criterion + dprime);
                  gradient = -( dnorm_at*w*(-n_responses_yes+n_responses*pnorm_at) );
                  gradient = gradient / ( pnorm_at*(-1+pnorm_at) );
              }
              break;

            case parameter_corr_mrsat:
                // No derivative for this parameter yet.
                break;


          };
          double x_Dtransform = coefs.TransformFn( coef_idx, CCoefs::FnConstrainDerivative );
          gradient = gradient * x_Dtransform;
          gradients[coef_idx] += gradient;
      }
    }
}


int CDesignMatrix::FindColumnIndex(std::string& column_name) {
    for(size_t i=0; i < mCoefNames.size(); i++){
        if(mCoefNames[i] == column_name)
          return (int)i;
    }
    return -1;
}

void CDesignMatrix::DetermineZeroRows(LogicalVector& zero_columns, std::vector<bool>& row_selected, bool all)
{
    for(size_t i=0; i < row_selected.size(); i++) {
      if(all) row_selected[i] = true;
      else    row_selected[i] = false;
    }
    
    CharacterVector column_names = zero_columns.names();
    for(int i=0; i < zero_columns.size(); i++) 
    {
      std::string column_name = as<std::string>(column_names[i]);
      bool zero_column = zero_columns[i];
      int col_idx = FindColumnIndex(column_name);

      for(int row_idx=0; row_idx < mDM.nrow(); row_idx++) 
      {
        double cur_value = mDM(row_idx, col_idx);
        if(all) {
            if(zero_column) row_selected[row_idx] = row_selected[row_idx] & (cur_value == 0.0);
            else            row_selected[row_idx] = row_selected[row_idx] & (cur_value != 0.0);
        } else {
            if(zero_column) row_selected[row_idx] = row_selected[row_idx] | (cur_value == 0.0);
            else            row_selected[row_idx] = row_selected[row_idx] | (cur_value != 0.0);
        }
      }
    }
}
    


/******************************
 *      class CCoefs          *
 ******************************/


CCoefs::CCoefs(DoubleVector& coefs, Function fn, bool use_names,
              CCoefConstraints* constraints) 
{
  mConstraints = constraints;

  if(fn == FnConstrain) {
    mUnconstrainedCoefs = coefs;
    mConstrainedCoefs = Constrain( coefs, use_names);

  } else {
    mConstrainedCoefs = coefs;
    mUnconstrainedCoefs = Unconstrain( coefs );
  }
}


double CCoefs::TransformWithOneBoundary(double x, double lower, Function fn) {
  switch(fn) {
    case FnConstrain: // [-Inf, +Inf] -> [lower, +Inf]
      return pow(x, 2)+lower;

    case FnConstrainDerivative:
      return 2*x;

    case FnUnconstrain:  // [lower, +Inf] -> [-Inf, +Inf]
      return sqrt(x - lower);
  };
  return nan("");
}

double CCoefs::TransformWithTwoBoundaries(double x, double lower, double upper, Function fn) {
  switch(fn) {
    case FnConstrain: // [-Inf, +Inf] -> [0, 1] -> [lower, upper]
    {
        double y = exp(x)/(1+exp(x));
        double z = y*(upper-lower)+lower;
        return(z);
    }
    break;

    case FnConstrainDerivative:
    {
      return exp(x)*(upper-lower)/(1+exp(x)) - pow(exp(x),2)*(upper-lower)/pow(1+exp(x), 2);
    }
    break;
    
    case FnUnconstrain: // [lower, upper] -> [0, 1] -> [-Inf, +Inf]
    {
        double y = (x-lower)/(upper-lower);
        double z = log(y/(1-y));
        return(z);
    }
    break;
  };
  return nan("");
}


double CCoefs::TransformFn(double raw_coef, double lower, double upper, Function fn)
{
  _dbg_method("TransformFn 1", 1);

  if( lower == upper ) {
    switch(fn) {
      case FnConstrainDerivative:
        return 0;

      case FnConstrain:
      case FnUnconstrain:
        return raw_coef;
    };
  
  } else if(lower == R_NegInf && upper == R_PosInf) {
    switch(fn) {
      case FnConstrainDerivative:
        return 1;

      case FnConstrain:
      case FnUnconstrain:
        return raw_coef;
    };

  } else if(upper == R_PosInf ) {
    return TransformWithOneBoundary(raw_coef, lower, fn);

  } else if(lower == R_NegInf) {
    return -1*TransformWithOneBoundary(-raw_coef, -upper, fn);
    
  } else {
    return TransformWithTwoBoundaries(raw_coef, lower, upper, fn);
    
  }
  return nan("");
}


double CCoefs::TransformFn(int coef_index, Function fn)
{
  _dbg_method("TransformFn 2", 1);
  switch(fn) {
    case FnConstrain:
    case FnConstrainDerivative:
        return TransformFn(mUnconstrainedCoefs[coef_index], mConstraints->mCoefsLower[coef_index], 
                         mConstraints->mCoefsUpper[coef_index], fn);

    case FnUnconstrain:
        return TransformFn(mConstrainedCoefs[coef_index], mConstraints->mCoefsLower[coef_index], 
                         mConstraints->mCoefsUpper[coef_index], fn);
  };
  return nan("");
}


Rcpp::DoubleVector CCoefs::Constrain(Rcpp::DoubleVector& raw_coefs, bool use_names) {
  _dbg_method("Constrain", 1);

  DoubleVector coefs; double val;

  for(int i=0; i < mConstraints->mCoefsLower.size(); i++)
  {
    val = TransformFn(raw_coefs[i], mConstraints->mCoefsLower[i], mConstraints->mCoefsUpper[i], FnConstrain);

    if(use_names) coefs.push_back( val, mConstraints->mCoefNames[i] );  
    else          coefs.push_back( val );  
  }
  return coefs;
}

DoubleVector CCoefs::Unconstrain(Rcpp::DoubleVector& raw_coefs)
{
  _dbg_method("Unconstrain", 1);

  DoubleVector coefs; double val;

  for(int i=0; i < raw_coefs.size(); i++) {
      val = TransformFn(raw_coefs[i], mConstraints->mCoefsLower[i], mConstraints->mCoefsUpper[i], FnUnconstrain);
      coefs.push_back( val );
  }
  return coefs;
}


double CCoefs::CoefsLL()
{
/*  _dbg_method("CoefsLL", 1);

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
  */
  return 0;
}




/******************************
 *   class CCoefConstraints   *
 ******************************/

CCoefConstraints::CCoefConstraints(Rcpp::NumericMatrix& coef_constraints) {
    _dbg_method("constructor", 1);
  
    mCoefsLower = coef_constraints(_, 0);
    mCoefsUpper = coef_constraints(_, 1);

    // make a copy for SetCoefValues() / ResetCoefValues()
    mCoefsLowerOriginal = clone( mCoefsLower );
    mCoefsUpperOriginal = clone( mCoefsUpper );

    List dimnames(NumericMatrix(coef_constraints).attr("dimnames"));
    if (dimnames.size() > 0) {
      RObject names = dimnames[0];
      if (!names.isNULL()) 
      {
        CharacterVector coef_names = CharacterVector(names);
        _dbg((0, "copying <%d> names", coef_names.size() ));
        for(int i=0; i < coef_names.size(); i++) {
          mCoefNames.push_back( as< std::string >(coef_names[i]) );
        }
      }
    }
}

CCoefConstraints::~CCoefConstraints() {
    _dbg_method("destructor", 1);
}

CCoefs CCoefConstraints::Constrain(Rcpp::DoubleVector& coefs, bool use_names) {
    return CCoefs(coefs, CCoefs::FnConstrain, use_names, this);
}

CCoefs CCoefConstraints::Unconstrain(Rcpp::DoubleVector& coefs) {
    return CCoefs(coefs, CCoefs::FnUnconstrain, false, this);
}

bool CCoefConstraints::SetCoefValue(std::string name, double value) {
    int coef_idx = FindCoefIndex(name);
    if(coef_idx == -1)
      return false;

//printf("setting %s to %.2f\n", name.c_str(), value);
    mCoefsLower[coef_idx] = value;
    mCoefsUpper[coef_idx] = value;
    return true;
}

bool CCoefConstraints::ResetCoefValue(std::string name) {
    int coef_idx = FindCoefIndex(name);
    if(coef_idx == -1)
      return false;

//printf("resetting %s to [%.2f, %.2f]\n", name.c_str(),  mCoefsLowerOriginal[coef_idx], mCoefsUpperOriginal[coef_idx]);
    mCoefsLower[coef_idx] = mCoefsLowerOriginal[coef_idx];
    mCoefsUpper[coef_idx] = mCoefsUpperOriginal[coef_idx];
    return true;
}

int CCoefConstraints::FindCoefIndex(std::string& coef_name) {
  for(size_t i=0; i < mCoefNames.size(); i++) {
    if(mCoefNames[i] == coef_name)
      return i;
  }
  return -1;
}

void CCoefConstraints::SetCoefValues(DoubleVector& values) {
    CharacterVector names = values.names();
    for(int i=0; i < values.size(); i++) {
        std::string coef_name = as<std::string>(names[i]);
        SetCoefValue(coef_name, values[i]);
    }
}

void CCoefConstraints::ResetCoefValues(CharacterVector& names) {
    for(int i=0; i < names.size(); i++) {
        std::string coef_name = as<std::string>(names[i]);
        ResetCoefValue(coef_name);
    }
}



/******************************
 *   class CSATFData          *
 ******************************/

CSATFData::CSATFData(Rcpp::CharacterVector& dv, Rcpp::NumericMatrix& dm,
              Rcpp::IntegerVector& dm_ncoef, Rcpp::NumericMatrix& constraints, 
              Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames):
                mDM(dm, dm_ncoef), mDV(dv, cnames, data), 
                mCoefConstraints(constraints)
{
    _dbg_method("constructor", 1);

    std::string name_time = as< std::string >(cnames["time"]);
    mTime = as< std::vector<double> >(data[name_time]);
    mEnabled = std::vector<bool>(mTime.size(), true);
}

CSATFData::~CSATFData() {
}



inline void save_LL(DoubleVector& LLVector, double cur_LL, bool by_row) {
    if(by_row) 
        LLVector.push_back( cur_LL );      
    else
        LLVector[0] += cur_LL;
}

inline bool valid_dprime(double dprime, DoubleVector& LLVector, bool by_row) {
    if( isnan(dprime) ) {
        printf("SATF WARNING: dprime is undefined.\n");
        save_LL(LLVector, R_NaN, by_row);
        return false;
    }
    return true;
}

inline bool valid_probability(double probability, DoubleVector& LLVector, bool by_row) {
    if(probability < 0 || probability > 1 || isnan(probability)) {
        printf("WARNING: approximation returned invalid probability: p_No = %f\n", probability);
        save_LL(LLVector, R_NaN, by_row);
        return false;
    }
    return true;
}


/* coefs: the coefs for computation
 * by_row: if true, returns a vector with log-likelihoods for every datapoint
 * tolerate_imprecisions: if true, tolerates and fixes minor imprecisions in the approximation to the conditional multivariate normal CDF
*/
DoubleVector CSATFData::ObjectiveFunction_Binary(CCoefs& coefs, bool by_row, bool tolerate_imprecisions, DoubleVector* gradient)
{
  _dbg_method("ObjectiveFunction_Binary", -1);
  assert( mDV.GetRVType() == rv_binary);
  _dbg((0, "rv_type: %d, datapoints: %d", mDV.GetRVType(), mDM.nDatapoints()));

  DoubleVector LLVector;
  double corr_mrsat = 0;

  if(!by_row)
    LLVector.push_back( 0 );

  // determine log-likelihood for each datapoint
  for(size_t idx = 0; idx < mDM.nDatapoints(); idx++)
  {
    if(!mEnabled[idx])
      continue;
    
    // get mrsat-correlation parameter if there is one
    if( mDM.HasParameter( CDesignMatrix::parameter_corr_mrsat ) ) {
      corr_mrsat = mDM.GetParameter(CDesignMatrix::parameter_corr_mrsat, coefs, idx);
      if( isnan(corr_mrsat) ) {
          printf("SATF WARNING: corr_mrsat is NaN.\n");
          save_LL(LLVector, R_NaN, by_row);
          return false;
      }
    }
      
    double criterion = mDM.ComputeCriterion(coefs, idx, mTime[idx]);
    double dprime = mDM.ComputeDprime(coefs, idx, mTime[idx]);

    if( !valid_dprime(dprime, LLVector, by_row) ) {
      mDM.ComputeCriterion(coefs, idx, mTime[idx], true);
      mDM.ComputeDprime(coefs, idx, mTime[idx], true);
      return LLVector;
    }
    
    _dbg((0, "idx <%d>, time <%.2f>, crit <%.2f>, dprime <%.2f>, resp. <%d>", 
              idx, mTime[idx], mDM.ComputeCriterion(coefs, idx, mTime[idx]), 
              mDM.ComputeDprime(coefs, idx, mTime[idx]), 
              mDV.ResponseYes(idx) ));
    
    int response = mDV.ResponseYes(idx);

    double cur_LL;
    if( corr_mrsat == 0.0 ) {
      cur_LL = mDM.ComputeLogLik_Independent(coefs, idx, mTime[idx], response, 1, gradient);

    } else 
    {
      double last_criterion = mDM.ComputeCriterion(coefs, idx-1, mTime[idx-1]);
      double last_dprime = mDM.ComputeDprime(coefs, idx-1, mTime[idx-1]);
      int last_response = mDV.ResponseYes(idx-1);

      // TODO: take care of last_dprime, and last_crit
      double crit_minus_psi = criterion - dprime;
      double last_crit_minus_psi = last_criterion - last_dprime;
      double pNo = pnorm_conditional(corr_mrsat, crit_minus_psi, last_crit_minus_psi, last_response==1, tolerate_imprecisions);

      _dbg((0, "corr.mrsat <%f>, crit_minus_psi <%.3f>, last_crit_minus_psi <%.3f>", corr_mrsat, crit_minus_psi, last_crit_minus_psi));

      if( !valid_probability(pNo, LLVector, by_row) ) {
        printf("coefs <%f, %f, %f, %f>\n", coefs[3], coefs[4], coefs[5], coefs[6]);
        mDM.ComputeCriterion(coefs, idx, mTime[idx], true);
        mDM.ComputeDprime(coefs, idx, mTime[idx], true);
        printf("<%f, %f, %f, %f>\n", criterion, dprime, last_criterion, last_dprime);
        printf("corr.mrsat <%f>, crit_minus_psi <%.3f>, last_crit_minus_psi <%.3f>\n", corr_mrsat, crit_minus_psi, last_crit_minus_psi);
        return LLVector;
      }
      
      if(mDV.ResponseYes(idx) == 1) cur_LL = log(1-pNo);
      else                          cur_LL = log(pNo);
    }
    save_LL(LLVector, cur_LL, by_row);
  }
  return LLVector;
}

DoubleVector CSATFData::ObjectiveFunction_Aggregate(CCoefs& coefs, bool by_row, DoubleVector* gradient)
{
  _dbg_method("ObjectiveFunction_Aggregate", 1);
  assert( mDV.GetRVType() == rv_aggregate);

  DoubleVector LLVector;
  if(!by_row)
    LLVector.push_back( 0 );

  for(size_t idx=0; idx < mDM.nDatapoints(); idx++)
  {
    if(!mEnabled[idx])
      continue;

    double cur_LL = mDM.ComputeLogLik_Independent(coefs, idx, mTime[idx], mDV.NResponsesYes(idx), mDV.NResponses(idx), gradient);

    _dbg((0, "idx <%d>, criterion <%.3f>, dprime <%.3f>, log-lik(%d/%d)=%.3f", 
              idx, mDM.ComputeCriterion(coefs, idx, mTime[idx]), 
              mDM.ComputeDprime(coefs, idx, mTime[idx]), mDV.NResponsesYes(idx), mDV.NResponses(idx), cur_LL));

    save_LL(LLVector, cur_LL, by_row);
  }
  _dbg((0, "returning <%f>", LLVector[0]));
  return LLVector;
}

Rcpp::DoubleVector CSATFData::ObjectiveFunction(DoubleVector& raw_coefs, bool by_row, 
                                                bool tolerate_imprecisions, DoubleVector* gradient)
{
  _dbg_method("ObjectiveFunction", 1);

  assert(coefs.length() == mnCoefs || coefs.length() == mnCoefs+1 );

  // transform coefficients
  CCoefs coefs = mCoefConstraints.Constrain(raw_coefs, false);

  // init the log-likelihood vector with the log-likelihood of the parameters, given the hyperparameters (or 0 if there are none)
  DoubleVector LL(1);
  LL[0] = coefs.CoefsLL();
  
  // compute the log-likelihood of the data
  switch(mDV.GetRVType()) 
  {
      case rv_aggregate:
        if(by_row) return      ObjectiveFunction_Aggregate(coefs, by_row, gradient);
        else       return LL + ObjectiveFunction_Aggregate(coefs, by_row, gradient);
        break;
        
      case rv_binary:
        if(by_row) return      ObjectiveFunction_Binary(coefs, by_row, tolerate_imprecisions, gradient);
        else       return LL + ObjectiveFunction_Binary(coefs, by_row, tolerate_imprecisions, gradient);
        break;
        
      case rv_dprime:
        // TODO: Not implemented for now.
        break;
  }
  return(LL);
}


void CSATFData::SelectSubset(LogicalVector& columns_zero, bool all) {
    mDM.DetermineZeroRows(columns_zero, mEnabled, all);
}

void CSATFData::ResetSubset(){
  for(size_t i=0; i < mEnabled.size(); i++) {
    mEnabled[i] = true;
  }
}
