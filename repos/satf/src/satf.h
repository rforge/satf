#ifndef __SATF_DATA_H__
#define __SATF_DATA_H__
#include <Rcpp.h>

//#define DEBUG
#include "debug.h"

enum RVType {
  rv_binary = 1,
  rv_aggregate = 2,
  rv_dprime = 3
};


class CResponseVariable {
  public:
    CResponseVariable(DoubleVector& predicted_criterion, 
                          CharacterVector& dv, CharacterVector& cnames,
                          DataFrame& data);
                          
    inline RVType GetRVType() {return mRVType; }
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
    RVType mRVType;
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

class CDesignMatrix
{
  public:
    typedef enum {
      parameter_asymptote = 0,
      parameter_invrate  = 1,
      parameter_intercept  = 2,
      parameter_invalid = 3
    } Parameter;
    
  public:
    CDesignMatrix(Rcpp::List& contrasts, Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames);
    inline size_t nCoefs() { return mDM.ncol(); }
    inline size_t nDatapoints() {return mDM.nrow(); }
    
    double ComputeDprime(Rcpp::DoubleVector& coefs, int datapoint_index, double time);

    double GetParameter(Parameter parameter, int datapoint_index, Rcpp::DoubleVector& coefs);
    
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

class CCoefConstraints {
  public:
    CCoefConstraints(Rcpp::NumericMatrix& coef_constraints);

    Rcpp::DoubleVector Constrain(Rcpp::DoubleVector& coefs, bool use_names);
    Rcpp::DoubleVector Unconstrain(Rcpp::DoubleVector& coefs);
    
    void AddCoefficient(double lower, double upper, std::string& name);

    double CoefsLL(Rcpp::DoubleVector& coefs);
    
  private:
    double TransformCoef(double raw_coefs, int i, bool constrain);

    double TransformWithOneBoundary(double x, double lower, bool constrain);
    double TransformWithTwoBoundaries(double x, double lower, double upper, bool constrain);
    
  private:
    Rcpp::DoubleVector mCoefsLower;
    Rcpp::DoubleVector mCoefsUpper;
    std::vector<std::string> mCoefNames;
    
    Rcpp::NumericMatrix mHyperparams;

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
    
    Rcpp::DoubleVector ObjectiveFunction(Rcpp::DoubleVector& coefs, bool by_row=false);
    
    Rcpp::DoubleVector ConstrainCoefs(Rcpp::DoubleVector& coefs, bool use_names);
    Rcpp::DoubleVector UnconstrainCoefs(Rcpp::DoubleVector& coefs);
    
    void AddCoefficient(double lower, double upper, std::string& name);

  private:
    Rcpp::DoubleVector ObjectiveFunction_Binary(Rcpp::DoubleVector& coefs, bool by_row=false);
    Rcpp::DoubleVector ObjectiveFunction_Aggregate(Rcpp::DoubleVector& coefs, bool by_row=false);

//  private:
  public:
    CDesignMatrix mDM;
    CResponseVariable mRV;
    CCoefConstraints mCoefConstraints;
    
    std::vector<double> mTime;
    
    _dbg_class_init;
};

inline Rcpp::DoubleVector CSATFData::ConstrainCoefs(Rcpp::DoubleVector& coefs, bool use_names) {
  return mCoefConstraints.Constrain(coefs, use_names);
}
inline Rcpp::DoubleVector CSATFData::UnconstrainCoefs(Rcpp::DoubleVector& coefs) {
  return mCoefConstraints.Unconstrain(coefs);
}
inline void CSATFData::AddCoefficient(double lower, double upper, std::string& name) {
  return mCoefConstraints.AddCoefficient(lower, upper, name);
}


DoubleVector ObjectiveFunctionCriterion(DoubleVector& coefs, DoubleVector& time, 
                                        IntegerVector& response, IntegerVector& trial_id);
DoubleVector ComputeCriterion(DoubleVector& coefs, DoubleVector& time);

#endif // __SATF_DATA_H__
