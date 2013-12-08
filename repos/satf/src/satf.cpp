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

#define GET(TYPE, NAME) Rcpp::clone( Rcpp::as<TYPE>(params[NAME]) )

#define debug_level 10

_dbg_class_set(CCoefConstraints, "CCoefConstraints", 0);
_dbg_class_set(CDependentVariable, "CDependentVariable", 0);
_dbg_class_set(CDesignMatrix, "CDesignMatrix",0);
_dbg_class_set(CSATFData, "CSATFData", 0);


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
  _dbg_function("constructor", 1);

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

/**********************************
 *   class CDesignMatrix  *
 **********************************/

CDesignMatrix::CDesignMatrix(Rcpp::NumericMatrix& dm, Rcpp::IntegerVector& dm_ncoef):  
                                        mCoefIndexFirst(parameter_invalid, -1),
                                        mCoefIndexLast(parameter_invalid, -1)
{
  _dbg_function("constructor", 1);

  mDM = Rcpp::clone(dm);

  CharacterVector ncoef_names = dm_ncoef.names();
  int cur_coef_index = 0;
  for(int i=0; i < dm_ncoef.size(); i++) {
    std::string cur_name = as<std::string>(ncoef_names[i]);
    Parameter parameter = StringToParameter( cur_name );
    mCoefIndexFirst[parameter] = cur_coef_index;
    mCoefIndexLast[parameter] = mCoefIndexFirst[parameter] + dm_ncoef[i] - 1;
    cur_coef_index += dm_ncoef[i];
  }
}

// void CDesignMatrix::StringToParameter(std::string& name) {  
// }


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

double CDesignMatrix::GetParameter(Parameter parameter, DoubleVector& coefs, int datapoint_index)  
{
    _dbg_function("GetParameter", 1);
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

double CDesignMatrix::ComputeDprime(Rcpp::DoubleVector& coefs, int datapoint_index, double time, bool log) 
{
  _dbg_function("ComputeDprime", 1);
  double asymptote = GetParameter(parameter_satf_asymptote, coefs, datapoint_index);
  if(asymptote == 0.0)
    return 0.0;
  double invrate   = GetParameter(parameter_satf_invrate,  coefs, datapoint_index);
  double intercept = GetParameter(parameter_satf_intercept, coefs, datapoint_index);
  if(log)
    printf("SATF(%f, %f, %f) = %f\n", asymptote, invrate, intercept, SATF(time, asymptote, invrate, intercept) );
  return SATF(time, asymptote, invrate, intercept);
}

// This is what Wickens (2001) calls lambda_center in equation (2.5) 
double CDesignMatrix::ComputeCriterion(Rcpp::DoubleVector& coefs, int datapoint_index, double time, bool log) 
{
  _dbg_function("ComputeCprime", 1);
  double min = GetParameter(parameter_bias_min, coefs, datapoint_index);
  double max = GetParameter(parameter_bias_max, coefs, datapoint_index);
  double invrate = GetParameter(parameter_bias_invrate, coefs, datapoint_index);
  double intercept = GetParameter(parameter_bias_intercept, coefs, datapoint_index);
  if(log)
    printf("criterion(%f, %f, %f, %f) = %f\n", min, max, invrate, intercept, (SATF(time, max-min, invrate, intercept) + min) );
  return SATF(time, max-min, invrate, intercept) + min;
}



/******************************
 *   class CCoefConstraints   *
 ******************************/

CCoefConstraints::CCoefConstraints(Rcpp::NumericMatrix& coef_constraints) {
    _dbg_function("constructor", 1);
  
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
    
  // FUTURE-TODO: Align the hyperparams by name
  // if( params.containsElementNamed("hyperparams") )
  //   mHyperparams = as<NumericMatrix>(params["hyperparams"]);
}

/*
void CCoefConstraints::AddCoefficient(double lower, double upper, std::string& name) {
    mCoefsLower.push_back(lower);
    mCoefsUpper.push_back(upper);
    mCoefNames.push_back(name);
}
*/

double CCoefConstraints::TransformWithOneBoundary(double x, double lower, bool constrain) {
  if(constrain) // [-Inf, +Inf] -> [lower, +Inf]
    return pow(x, 2)+lower;
  else // [lower, +Inf] -> [-Inf, +Inf]
    return sqrt(x - lower);
}

double CCoefConstraints::TransformWithTwoBoundaries(double x, double lower, double upper, bool constrain) {
  if(constrain) { // [-Inf, +Inf] -> [0, 1] -> [lower, upper]
    double y = exp(x)/(1+exp(x));
    double z = y*(upper-lower)+lower;
    return(z);

  } else // [lower, upper] -> [0, 1] -> [-Inf, +Inf]
  { 
    double y = (x-lower)/(upper-lower);
    double z = log(y/(1-y));
    return(z);
  }
}
double CCoefConstraints::TransformCoef(double raw_coef, int i, bool constrain)
{
  _dbg_function("TransformCoef", 1);

  if( IsCoefFixed(i) ) {
    return raw_coef;

  } else if(mCoefsLower[i] == R_NegInf && mCoefsUpper[i] == R_PosInf) {
    return raw_coef; 

  } else if(mCoefsUpper[i] == R_PosInf ) {
    return TransformWithOneBoundary(raw_coef, mCoefsLower[i], constrain);

  } else if( mCoefsLower[i] == R_NegInf) {
    return -1*TransformWithOneBoundary(-raw_coef, -mCoefsUpper[i], constrain);
    
  } else {
    return TransformWithTwoBoundaries(raw_coef, mCoefsLower[i], mCoefsUpper[i], constrain);
  }
}

Rcpp::DoubleVector CCoefConstraints::Constrain(Rcpp::DoubleVector& raw_coefs, bool use_names) {
  _dbg_function("Constrain", 1);

  DoubleVector coefs;
  double val;

  for(int i=0; i < mCoefsLower.size(); i++)
  {
    val = TransformCoef(raw_coefs[i], i, true);

    if(use_names)
      coefs.push_back( val, mCoefNames[i] );  
    else
      coefs.push_back( val );  
  }
  return coefs;
}

DoubleVector CCoefConstraints::Unconstrain(Rcpp::DoubleVector& raw_coefs)
{
  _dbg_function("Unconstrain", 1);

  DoubleVector coefs;
  double val;

  for(int i=0; i < raw_coefs.size(); i++) {
      val = TransformCoef(raw_coefs[i], i, false);
      coefs.push_back( val );
  }
  return coefs;
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

double CCoefConstraints::CoefsLL(DoubleVector& coefs)
{
/*  _dbg_function("CoefsLL", 1);

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
 *   class CSATFData          *
 ******************************/

CSATFData::CSATFData(Rcpp::CharacterVector& dv, Rcpp::NumericMatrix& dm,
              Rcpp::IntegerVector& dm_ncoef, Rcpp::NumericMatrix& constraints, 
              Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames):
                mDM(dm, dm_ncoef), mDV(dv, cnames, data), 
                mCoefConstraints(constraints)
{
  _dbg_function("constructor", 1);

    std::string name_time = as< std::string >(cnames["time"]);
    mTime = as< std::vector<double> >(data[name_time]);
    mDisabled = std::vector<bool>(mTime.size(), false);
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
DoubleVector CSATFData::ObjectiveFunction_Binary(DoubleVector& coefs, bool by_row, bool tolerate_imprecisions)
{
  _dbg_function("ObjectiveFunction_Binary", -1);
  assert( mDV.GetRVType() == rv_binary);
  _dbg((0, "rv_type: %d, datapoints: %d", mDV.GetRVType(), mDM.nDatapoints()));

  DoubleVector LLVector;
  double corr_mrsat = 0;
  double last_dprime = 0, last_criterion = 0;
  int last_response = 0;

  if(!by_row)
    LLVector.push_back( 0 );

  // determine log-likelihood for each datapoint
  for(size_t idx = 0; idx < mDM.nDatapoints(); idx++)
  {
    if(mDisabled[idx])
      continue;

    int response = mDV.ResponseYes(idx);
    double dprime = mDM.ComputeDprime(coefs, idx, mTime[idx]);
    double criterion = mDM.ComputeCriterion(coefs, idx, mTime[idx]);

    if( !valid_dprime(dprime, LLVector, by_row) ) {
      printf("criterion %f\n", mDM.ComputeCriterion(coefs, idx, mTime[idx], true) );
      printf("dprime %f\n", mDM.ComputeDprime(coefs, idx, mTime[idx], true) );
      return LLVector;
    }
    
    _dbg((0, "idx <%d>, time <%.2f>, crit <%.2f>, dprime <%.2f>, resp. <%d>", idx, mTime[idx], criterion, dprime, response ));
    
    // get mrsat-correlation parameter if there is one
    if( mDM.HasParameter( CDesignMatrix::parameter_corr_mrsat ) ) {
      corr_mrsat = mDM.GetParameter(CDesignMatrix::parameter_corr_mrsat, coefs, idx);
      if( isnan(dprime) ) {
          printf("SATF WARNING: corr_mrsat is NaN.\n");
          save_LL(LLVector, R_NaN, by_row);
          return false;
      }
    }

    double pNo, cur_LL;
    if( corr_mrsat == 0.0 ) {
        cur_LL = _pnorm(criterion, dprime, 1.0, response==0, true);

    } else 
    {
      double crit_minus_psi = criterion - dprime;
      double last_crit_minus_psi = last_criterion - last_dprime;
      pNo = pnorm_conditional(corr_mrsat, crit_minus_psi, last_crit_minus_psi, last_response==1, tolerate_imprecisions);

      _dbg((0, "corr.mrsat <%f>, crit_minus_psi <%.3f>, last_crit_minus_psi <%.3f>", corr_mrsat, crit_minus_psi, last_crit_minus_psi));

      if( !valid_probability(pNo, LLVector, by_row) ) {
        printf("coefs <%f, %f, %f, %f>\n", coefs[3], coefs[4], coefs[5], coefs[6]);
        printf("criterion %f\n", mDM.ComputeCriterion(coefs, idx, mTime[idx], true) );
        printf("dprime %f\n", mDM.ComputeDprime(coefs, idx, mTime[idx], true) );
        printf("<%f, %f, %f, %f>\n", criterion, dprime, last_criterion, last_dprime);
        printf("corr.mrsat <%f>, crit_minus_psi <%.3f>, last_crit_minus_psi <%.3f>\n", corr_mrsat, crit_minus_psi, last_crit_minus_psi);
        return LLVector;
      }
      
      if(mDV.ResponseYes(idx) == 1) cur_LL = log(1-pNo);
      else                          cur_LL = log(pNo);
    }
    last_dprime = dprime;
    last_criterion = criterion;
    last_response = response;

    save_LL(LLVector, cur_LL, by_row);
  }
  return LLVector;
}

DoubleVector CSATFData::ObjectiveFunction_Aggregate(Rcpp::DoubleVector& coefs, bool by_row)
{
  _dbg_function("ObjectiveFunction_Aggregate", -1);
//printf("coef 2 <%f>\n", coefs[1]);

  assert( mDV.GetRVType() == rv_aggregate);

  DoubleVector LLVector;
  if(!by_row)
    LLVector.push_back( 0 );

  for(size_t idx=0; idx < mDM.nDatapoints(); idx++)
  {
    if(mDisabled[idx])
      continue;

    double dprime = mDM.ComputeDprime(coefs, idx, mTime[idx]);
    double criterion = mDM.ComputeCriterion(coefs, idx, mTime[idx]);

    double pYes = 1-_pnorm(criterion, dprime);
    int n_responses_yes = mDV.NResponsesYes(idx);
    int n_responses = mDV.NResponses(idx);
    double cur_LL = _dbinom(n_responses_yes, n_responses, pYes, true );
    _dbg((0, "idx <%d>, criterion <%.3f>, dprime <%.3f>, p_yes <%.3f>, log-lik(%d/%d)=%.3f", 
              idx, criterion, dprime, pYes, n_responses_yes, n_responses, cur_LL));

    save_LL(LLVector, cur_LL, by_row);
  }
  _dbg((0, "returning <%f>", LLVector[0]));
  return LLVector;
}

Rcpp::DoubleVector CSATFData::ObjectiveFunction(DoubleVector& raw_coefs, bool by_row, bool tolerate_imprecisions)
{
  _dbg_function("ObjectiveFunction", 1);

  assert(coefs.length() == mnCoefs || coefs.length() == mnCoefs+1 );

  // transform coefficients
  DoubleVector coefs = ConstrainCoefs(raw_coefs, false);

  // init the log-likelihood vector with the log-likelihood of the parameters, given the hyperparameters (or 0 if there are none)
  DoubleVector LL(1);
  LL[0] = mCoefConstraints.CoefsLL(coefs);
  
  // compute the log-likelihood of the data
  switch(mDV.GetRVType()) 
  {
      case rv_aggregate:
        if(by_row)
          return ObjectiveFunction_Aggregate(coefs, by_row);
        else
          return LL + ObjectiveFunction_Aggregate(coefs, by_row);
        break;
        
      case rv_binary:
        if(by_row)
          return ObjectiveFunction_Binary(coefs, by_row, tolerate_imprecisions);
        else
          return LL + ObjectiveFunction_Binary(coefs, by_row, tolerate_imprecisions);
        break;
        
      case rv_dprime:
        // TODO: Not implemented for now.
        break;
  }
  return(LL);
}


void CSATFData::SelectSubset(LogicalVector& selection) {
  for(size_t i=0; i < mDisabled.size(); i++) {
    mDisabled[i] = !selection[i];
  }
}
void CSATFData::ResetSubset(){
  for(size_t i=0; i < mDisabled.size(); i++) {
    mDisabled[i] = false;
  }
}