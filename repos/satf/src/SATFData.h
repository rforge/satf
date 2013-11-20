#ifndef __SATF_DATA_H__
#define __SATF_DATA_H__
#include <Rcpp.h>
#include "debug.h"

enum enumRVType {
  rv_binary = 1,
  rv_aggregate = 2,
  rv_dprime = 3
};


class CSATFResponseVariable {
  public:
    CSATFResponseVariable(DoubleVector& predicted_criterion, 
                          CharacterVector& dv, CharacterVector& cnames,
                          DataFrame& data);
                          
    inline enumRVType GetRVType() {return mRVType; }
    inline double PredictedCriterion(int idx) {return mPredictedCriterion[idx]; }

    // response variable 1
    inline int ResponseYes(int idx) {return mResponseYes[idx]; }
    inline int TrialId(int idx) {return mTrialId[idx]; }
    inline bool HasTrialId() {return mTrialId.size() > 0; }

    // response variable 2
    inline int ResponseDprime(int idx) {return mResponseDprime[idx]; }

    // response variable 3
    inline int NResponsesYes(int idx) {return mnResponsesYes[idx]; }
    inline int NResponses(int idx) {return mnResponses[idx]; }

//  private:
    public:
    enumRVType mRVType;
    Rcpp::DoubleVector mPredictedCriterion;

    // response variable 1
    std::vector<int> mResponseYes;
    std::vector<int> mTrialId;
    
    // response variable 2
    std::vector<double> mResponseDprime;
    
    // response variable 3
    IntegerVector mnResponsesYes;
    IntegerVector mnResponses;

    _dbg_class_init;
};

class CSATFDesignMatrix
{
  public:
    typedef enum {
      parameter_lambda = 0,
      parameter_beta  = 1,
      parameter_delta  = 2,
      parameter_invalid = 3
    } Parameter;
    
  public:
    CSATFDesignMatrix(Rcpp::List& contrasts, Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames);
    inline size_t nCoefs() { return mDM.ncol(); }
    inline size_t nDatapoints() {return mDM.nrow(); }
    
    double ComputeDprime(std::vector<double>& coefs, int datapoint_index, double time);

    double GetParameter(Parameter parameter, int datapoint_index, std::vector<double>& coefs);
    
    Parameter StringToParameter(std::string& name);

//  private:
    public:
    Rcpp::NumericMatrix mDM;
    
    std::vector<int> mCoefIndexFirst;
    std::vector<int> mCoefIndexLast;
    
    unsigned int mnCoefs;
    unsigned int mnDatapoints;

    _dbg_class_init;
};

class CSATFData 
{
  public:
    CSATFData(Rcpp::CharacterVector& dv, Rcpp::List& contrasts,
              Rcpp::NumericMatrix& coef_constraints, DataFrame& data, 
              Rcpp::DoubleVector& predicted_criterion, 
              Rcpp::CharacterVector& cnames);
    ~CSATFData();
    
    double ObjectiveFunction(Rcpp::DoubleVector& coefs);
    
    std::vector<double> TransformCoefs(Rcpp::DoubleVector& coefs, bool back=false);
    
  private:
    double TransformCoef(double raw_coefs, int i, bool back);
    double CoefsLL(std::vector<double>& coefs);

    double ObjectiveFunction_Binary(std::vector<double>& coefs);
    double ObjectiveFunction_Aggregate(std::vector<double>& coefs);

//  private:
  public:
    CSATFDesignMatrix mDM;
    CSATFResponseVariable mRV;
    std::vector<double> mTime;
    Rcpp::NumericMatrix mHyperparams;
    
    Rcpp::DoubleVector mParametersLower;
    Rcpp::DoubleVector mParametersUpper;
    
    std::vector<int> mDbgTime;
    
    _dbg_class_init;
};

#endif // __SATF_DATA_H__
