#ifndef __SATF_DATA_H__
#define __SATF_DATA_H__
#include <Rcpp.h>

enum enumRVType {
  rv_binary = 1,
  rv_aggregate = 2,
  rv_dprime = 3
};

class CSATFResponseVariable {

  public:
    CSATFResponseVariable(Rcpp::List& params);
    inline enumRVType GetRVType() {return mRVType; }
    inline double PredictedCriterion(int idx) {return mPredictedCriterion[idx]; }

    // response variable 1
    inline int ResponseYes(int idx) {return mResponseYes[idx]; }
    inline int TrialId(int idx) {return mTrialId[idx]; }
    inline bool HasTrialId() {return mTrialId.length() != 0; }

    // response variable 2
    inline int ResponseDprime(int idx) {return mResponseDprime[idx]; }

    // response variable 3
    inline int ResponseHits(int idx) {return mResponseHits[idx]; }
    inline int NGrammatical(int idx) {return mnGrammatical[idx]; }
    inline int ResponseFAs(int idx) {return mResponseFAs[idx]; }
    inline int NUngrammatical(int idx) {return mnUngrammatical[idx]; }

  private:
    enumRVType mRVType;
    Rcpp::DoubleVector mPredictedCriterion;

    // response variable 1
    Rcpp::IntegerVector mResponseYes;
    Rcpp::IntegerVector mTrialId;
    
    // response variable 2
    Rcpp::DoubleVector mResponseDprime;
    
    // response variable 3
    Rcpp::IntegerVector mResponseHits;
    Rcpp::IntegerVector mnGrammatical;
    Rcpp::IntegerVector mResponseFAs;
    Rcpp::IntegerVector mnUngrammatical;
};

class CSATFDesignMatrix 
{
  public:
    CSATFDesignMatrix(Rcpp::List& params);
    inline int NCoefs() { return mDM.nrow(); }
    inline int NDatapoints() {return mDM.ncol(); }
    
    double ComputeDprime(Rcpp::DoubleVector& coefs, int datapoint_index, double time);

    double GetParam(Rcpp::IntegerVector& rows, int datapoint_index, Rcpp::DoubleVector& coefs);  
  
  private:
    Rcpp::NumericMatrix mDM;
    
    Rcpp::IntegerVector mRowsLambda; 
    Rcpp::IntegerVector mRowsBeta; 
    Rcpp::IntegerVector mRowsDelta;
};

class CSATFData 
{
  public:
    CSATFData( Rcpp::List& params);
    double ObjectiveFunction(Rcpp::DoubleVector& coefs, Rcpp::DoubleVector& fixed_coefs);
    
  private:  
    void CoefsTransform(Rcpp::DoubleVector& coefs, Rcpp::IntegerVector& transform_vec, 
      		              Rcpp::DoubleVector& minimum, Rcpp::DoubleVector& maximum);
    double CoefsLL(Rcpp::DoubleVector& coefs);

    double ObjectiveFunction_Binary(DoubleVector& coefs);
    double ObjectiveFunction_Aggregate(DoubleVector& coefs);

  private:
    Rcpp::DoubleVector mFixedparams;
    
    CSATFDesignMatrix mDM;
    int mnCoefs;  int mnDatapoints;
    
    Rcpp::IntegerVector mTransformVec;
    Rcpp::DoubleVector mParamsMinimum;
    Rcpp::DoubleVector mParamsMaximum;
    
    Rcpp::DoubleVector mTime;

    Rcpp::NumericMatrix mHyperparams;

    CSATFResponseVariable mRV;
};

#endif // __SATF_DATA_H__
