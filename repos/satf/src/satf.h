#ifndef __SATF_DATA_H__
#define __SATF_DATA_H__
#include <Rcpp.h>

#include "debug.h"
#include "satf_math.h"


class CCoefs;

enum RVType {
  rv_binary = 1,
  rv_aggregate = 2,
  rv_dprime = 3
};

typedef enum EParameter {
  parameter_satf_asymptote = 0,
  parameter_satf_invrate  = 1,
  parameter_satf_intercept  = 2,
  parameter_bias_min = 3,
  parameter_bias_max  = 4,
  parameter_bias_invrate  = 5,
  parameter_bias_intercept  = 6,
  parameter_corr_mrsat  = 7,
  parameter_invalid = 8
} Parameter;

class CDesignMatrix
{
  public:
    CDesignMatrix(Rcpp::NumericMatrix& dm, Rcpp::IntegerVector& dm_ncoef);
    ~CDesignMatrix();
    
    inline size_t nCoefs() { return mDM.ncol(); }
    inline size_t nDatapoints() {return mDM.nrow(); }
    
    void DetermineZeroRows(Rcpp::LogicalVector& columns_zero, std::vector<bool>& row_selected, bool all);
    
    SNAEPoint ComputeDprime(CCoefs& coefs, int datapoint_index, double time, bool log=false);
    SNAEPoint ComputeCriterion(CCoefs& coefs, int datapoint_index, double time, bool log=false);

    double GetParameter(Parameter parameter, CCoefs& coefs, int datapoint_index);
    bool HasParameter(Parameter parameter);
    // std::vector<Parameter> GetParameterTypes();

    Parameter StringToParameter(std::string& name);

    inline double operator()(int i, int j);
    inline std::vector<int>& coef_types();
    
  private:
    int FindColumnIndex(std::string& column_name);

//  private:
    public:
    std::vector<std::string> mCoefNames;
    std::vector<int> mCoefTypes;

    std::vector<int> mCoefIndexFirst;
    std::vector<int> mCoefIndexLast;

    Rcpp::NumericMatrix mDM;

    unsigned int mnCoefs;
    unsigned int mnDatapoints;

    _dbg_class_init;
};

inline bool CDesignMatrix::HasParameter(Parameter parameter) {
    if(mCoefIndexFirst[parameter] == -1)
        return false;
    return true;
}

inline double CDesignMatrix::operator()(int i, int j) {
  return mDM(i, j);
}

inline std::vector<int>& CDesignMatrix::coef_types() {
  return mCoefTypes;
}


class CDataPoint {
  public:
    CDataPoint(int _n_responses_yes, int _n_responses, double _time, int _index, CDesignMatrix* dm);
  
    double ComputeLogLik(CCoefs& coefs);
    double ComputeLogLik(CCoefs& coefs, double corr_mrsat, CDataPoint& last_datapoint, bool tolerate_imprecisions=false);

    bool Update(CCoefs& coefs);
    
    void ResetDprime();
    void ResetCriterion();

    Rcpp::DoubleVector ComputeLogLikGradient(CCoefs& coefs);
    
private:
    bool UpdateDprime(CCoefs& coefs);
    bool UpdateCriterion(CCoefs& coefs);

public:
  SNAEPoint dprime;
  SNAEPoint criterion;
  int n_responses_yes;
  int n_responses;  
  double time;
  int index;  
  CDesignMatrix* mDM;
};

inline CDataPoint::CDataPoint(int _n_responses_yes, int _n_responses, double _time, int _index, CDesignMatrix* dm) {
   n_responses_yes = _n_responses_yes;
   n_responses = _n_responses;
   time = _time;
   index = _index;
   mDM = dm;
}


inline bool CDataPoint::Update(CCoefs& coefs) {
  if(!UpdateDprime(coefs))
    return false;
  if(!UpdateCriterion(coefs))
    return false;
  return true;
}

inline bool CDataPoint::UpdateDprime(CCoefs& coefs) {
  dprime = mDM->ComputeDprime(coefs, index, time);
  return !isnan(dprime.value);
}
inline bool CDataPoint::UpdateCriterion(CCoefs& coefs) {
  criterion = mDM->ComputeCriterion(coefs, index, time);
  return !isnan(criterion.value);
}

inline void CDataPoint::ResetDprime() {
  dprime.reset();
}
inline void CDataPoint::ResetCriterion() {
  criterion.reset();
}


class CCoefConstraints
{
  friend CCoefs;

  public:
    CCoefConstraints(Rcpp::NumericMatrix& coef_constraints);
    ~CCoefConstraints();
    
    CCoefs Constrain(Rcpp::DoubleVector& coefs, bool use_names);
    CCoefs Unconstrain(Rcpp::DoubleVector& coefs);

    void SetCoefValues(Rcpp::DoubleVector& values);
    void ResetCoefRanges(Rcpp::CharacterVector& names);

    bool IsCoefFixed(int coef_index);
    int  FindCoefIndex(std::string& name);
    
  private:
    bool SetCoefValue(std::string name, double value);
    bool ResetCoefRange(std::string name);

  protected:
  //private:
      Rcpp::DoubleVector mCoefsLower;
      Rcpp::DoubleVector mCoefsUpper;
      std::vector<std::string> mCoefNames;

      Rcpp::DoubleVector mCoefsLowerOriginal;
      Rcpp::DoubleVector mCoefsUpperOriginal;
      
      Rcpp::NumericMatrix mHyperparams;

    _dbg_class_init;
};

inline bool CCoefConstraints::IsCoefFixed(int coef_index) {
    return (mCoefsLower[coef_index] == mCoefsUpper[coef_index]);
}


class CCoefs {
  public:
    typedef enum EFunction {
      FnConstrain,
      FnConstrainDerivative,
      FnUnconstrain
    } Function;

public:
    CCoefs(Rcpp::DoubleVector& unconstrained_coefs, Function fn, 
           bool use_names, CCoefConstraints* constraints);

public:

    double CoefsLL();

    Rcpp::DoubleVector constrained();
    Rcpp::DoubleVector unconstrained();

    Rcpp::DoubleVector Constrain(Rcpp::DoubleVector& coefs, bool use_names);
    Rcpp::DoubleVector Unconstrain(Rcpp::DoubleVector& coefs);

    double TransformFn(int coef_index, Function fn);

    inline int size() { return mConstrainedCoefs.size(); } 
    inline double operator[](int i) { return mConstrainedCoefs[i]; } 

private:
    double TransformFn(double x, double lower, double upper, Function fn);
    double TransformWithOneBoundary(double x, double lower, Function fn);
    double TransformWithTwoBoundaries(double x, double lower, double upper, Function fn);

private:
    Rcpp::DoubleVector mUnconstrainedCoefs;
    Rcpp::DoubleVector mConstrainedCoefs;
  
    CCoefConstraints* mConstraints;

    _dbg_class_init;
};

inline Rcpp::DoubleVector CCoefs::constrained() {
  return mConstrainedCoefs;
}

inline Rcpp::DoubleVector CCoefs::unconstrained() {
  return mUnconstrainedCoefs;
}



class CDataContainer
{
  public:
    CDataContainer(Rcpp::CharacterVector& dv, Rcpp::NumericMatrix& dm,
              Rcpp::IntegerVector& dm_ncoef, Rcpp::NumericMatrix& constraints, 
              Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames);
    ~CDataContainer();
    
    Rcpp::DoubleVector ObjectiveFunction(Rcpp::DoubleVector& coefs, bool by_row=false, 
                                         bool tolerate_imprecisions=false);
    Rcpp::DoubleVector ObjectiveFunctionGradient(Rcpp::DoubleVector& coefs,
                                         bool by_row=false, bool tolerate_imprecisions=false);
    
    Rcpp::DoubleVector ConstrainCoefs(Rcpp::DoubleVector& coefs, bool use_names);
    Rcpp::DoubleVector UnconstrainCoefs(Rcpp::DoubleVector& coefs);

    void SelectSubset(Rcpp::LogicalVector& zero_columns, bool all);
    void ResetSubset();

    void SetCoefValues(Rcpp::DoubleVector& values);
    void ResetCoefRanges(Rcpp::CharacterVector& names);

//  private:
  public:
    CDesignMatrix mDM;
    CCoefConstraints mCoefConstraints;
    std::vector<CDataPoint> mDatapoints;
    std::vector<bool> mEnabled;
    
    _dbg_class_init;
};

inline void CDataContainer::SetCoefValues(Rcpp::DoubleVector& values) {
   mCoefConstraints.SetCoefValues(values);
}

inline void CDataContainer::ResetCoefRanges(Rcpp::CharacterVector& names) {
   mCoefConstraints.ResetCoefRanges(names);
}

inline Rcpp::DoubleVector CDataContainer::ConstrainCoefs(Rcpp::DoubleVector& coefs, bool use_names) {
  return mCoefConstraints.Constrain(coefs, use_names).constrained();
}
inline Rcpp::DoubleVector CDataContainer::UnconstrainCoefs(Rcpp::DoubleVector& coefs) {
  return mCoefConstraints.Unconstrain(coefs).unconstrained();
}


#endif // __SATF_DATA_H__
