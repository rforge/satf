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


class CDependentVariable
{
    public:
      CDependentVariable(Rcpp::CharacterVector& dv, Rcpp::CharacterVector& cnames, Rcpp::DataFrame& data);
                            
      inline RVType GetRVType() {return mDVType; }

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
        RVType mDVType;

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
    
  public:
    CDesignMatrix(Rcpp::NumericMatrix& dm, Rcpp::IntegerVector& dm_ncoef);

    inline size_t nCoefs() { return mDM.ncol(); }
    inline size_t nDatapoints() {return mDM.nrow(); }
    
    double ComputeDprime(Rcpp::DoubleVector& coefs, int datapoint_index, double time, bool log=false);
    double ComputeCriterion(Rcpp::DoubleVector& coefs, int datapoint_index, double time, bool log=false);

    double GetParameter(Parameter parameter, Rcpp::DoubleVector& coefs, int datapoint_index);
    bool HasParameter(Parameter parameter);

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

inline bool CDesignMatrix::HasParameter(Parameter parameter) {
    if(mCoefIndexFirst[parameter] == -1)
        return false;
    return true;
}


class CCoefConstraints
{
  public:
    CCoefConstraints(Rcpp::NumericMatrix& coef_constraints);

    Rcpp::DoubleVector Constrain(Rcpp::DoubleVector& coefs, bool use_names);
    Rcpp::DoubleVector Unconstrain(Rcpp::DoubleVector& coefs);
    
    void SetCoefValues(DoubleVector& values);
    void ResetCoefValues(CharacterVector& names);

    double CoefsLL(Rcpp::DoubleVector& coefs);

    bool IsCoefFixed(int coef_index);
    int FindCoefIndex(std::string& name);
    
  private:
    bool SetCoefValue(std::string name, double value);
    bool ResetCoefValue(std::string name);

    double TransformCoef(double raw_coefs, int i, bool constrain);
    double TransformWithOneBoundary(double x, double lower, bool constrain);
    double TransformWithTwoBoundaries(double x, double lower, double upper, bool constrain);
    
  private:
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



class CSATFData 
{
  public:
    CSATFData(Rcpp::CharacterVector& dv, Rcpp::NumericMatrix& dm,
              Rcpp::IntegerVector& dm_ncoef, Rcpp::NumericMatrix& constraints, 
              Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames);
    ~CSATFData();
    
    Rcpp::DoubleVector ObjectiveFunction(Rcpp::DoubleVector& coefs, bool by_row=false, bool tolerate_imprecisions=false);
    
    Rcpp::DoubleVector ConstrainCoefs(Rcpp::DoubleVector& coefs, bool use_names);
    Rcpp::DoubleVector UnconstrainCoefs(Rcpp::DoubleVector& coefs);

    void SelectSubset(LogicalVector& selection);
    void ResetSubset();

    void SetCoefValues(DoubleVector& values);
    void ResetCoefValues(CharacterVector& names);


  private:
    Rcpp::DoubleVector ObjectiveFunction_Binary(Rcpp::DoubleVector& coefs, bool by_row=false, bool tolerate_imprecisions=false);
    Rcpp::DoubleVector ObjectiveFunction_Aggregate(Rcpp::DoubleVector& coefs, bool by_row=false);

//  private:
  public:
        CDesignMatrix mDM;
        CDependentVariable mDV;
        CCoefConstraints mCoefConstraints;
        
        std::vector<double> mTime;
        std::vector<bool> mDisabled;
    
    _dbg_class_init;
};

inline void CSATFData::SetCoefValues(DoubleVector& values) {
   mCoefConstraints.SetCoefValues(values);
}

inline void CSATFData::ResetCoefValues(CharacterVector& names) {
   mCoefConstraints.ResetCoefValues(names);
}

inline Rcpp::DoubleVector CSATFData::ConstrainCoefs(Rcpp::DoubleVector& coefs, bool use_names) {
  return mCoefConstraints.Constrain(coefs, use_names);
}
inline Rcpp::DoubleVector CSATFData::UnconstrainCoefs(Rcpp::DoubleVector& coefs) {
  return mCoefConstraints.Unconstrain(coefs);
}


#endif // __SATF_DATA_H__
