
#include <Rcpp.h>
#include "satf.h"
using namespace Rcpp;


/******************************
 *   class CDataContainer     *
 ******************************/

CDataContainer::CDataContainer(Rcpp::CharacterVector& dv, Rcpp::NumericMatrix& dm,
                  Rcpp::IntegerVector& dm_ncoef, Rcpp::NumericMatrix& constraints, 
                  Rcpp::DataFrame& data, Rcpp::CharacterVector& cnames):
                        mCoefConstraints(constraints, dm_ncoef)
{
    _dbg_method("constructor", 1);

    mObjectiveFunctionCalls = 0;

    // save unique rows of the design matrix with lists of
    // corresponding datapoint indices
    for(int row_idx=0; row_idx < dm.nrow(); row_idx++) {
      DoubleVector dm_row = dm(row_idx, _);
      _dbg((0, "LENGTH: %d", dm_row.size() ));
      AddDMRow(dm_row, row_idx);
    }
      
    std::string name_time = as< std::string >(cnames["time"]);
    DoubleVector time = as< DoubleVector >(data[ name_time ]);
  
    if( dv.containsElementNamed("response") ) {
      
        std::string name_response = as< std::string >(dv["response"]);
        IntegerVector response_yes = as< IntegerVector >(data[ name_response ]);
        CDataPoint* last_datapoint = NULL;
        for(int i=0; i < data.nrows(); i++) {
          CDataPoint* dpoint = new CDataPoint(response_yes[i], 1, time[i], i, last_datapoint);
          mDatapoints.push_back( dpoint );
          last_datapoint = dpoint;
        }
  
    } else if( dv.containsElementNamed("dprime") ) {
        // not implemented
  
    } else if( dv.containsElementNamed("n.responses") ) {
        std::string name_yes = as< std::string >(dv["n.responses.yes"]);
        std::string name_all = as< std::string >(dv["n.responses"]);
        IntegerVector n_responses_yes = as< IntegerVector >(data[ name_yes ]);
        IntegerVector n_responses = as< IntegerVector >(data[ name_all ]);
        CDataPoint* last_datapoint = NULL;
        for(int i=0; i < data.nrows(); i++) {
          CDataPoint* dpoint = new CDataPoint( n_responses_yes[i], n_responses[i], time[i], i, last_datapoint);
          mDatapoints.push_back( dpoint );
          last_datapoint = dpoint;
        }
    }
      
    mEnabled = std::vector<bool>(mUniqueDMRows.size(), true);
    
    mForceUpdate = true;
}

CDataContainer::~CDataContainer()
{
  for(size_t i=0; i < mDatapoints.size(); i++) {
      delete mDatapoints[i];
      mDatapoints[i] = NULL;
  }
}


int CDataContainer::FindDMRow(CDesignMatrixRow& row)
{
  _dbg_method("FindDMRow", 1);

  for(size_t j=0; j < mUniqueDMRows.size(); j++) {
    if(mUniqueDMRows[j] == row)
      return j;
  }
  return -1;
}

void CDataContainer::AddDMRow(DoubleVector& dm_row, int datapoint_index)
{
  _dbg_method("AddDMRow", 1);

  CDesignMatrixRow row(dm_row, mCoefConstraints);
  int idx = FindDMRow(row);
  if(idx == -1) {
    mUniqueDMRows.push_back(row);
    idx = mUniqueDMRows.size()-1;
  }
  mUniqueDMRows[idx].AddDatapointIndex(datapoint_index);
}


inline void save_LL(DoubleVector& LLVector, double cur_LL, int idx, bool by_row) {
    if(by_row) 
        LLVector[idx] = cur_LL;
    else
        LLVector[0] += cur_LL;
}

Rcpp::DoubleVector CDataContainer::ObjectiveFunction(DoubleVector& raw_coefs, bool by_row, 
                                                     bool tolerate_imprecisions, bool force_update)
{
  _dbg_method("ObjectiveFunction", 1);

  if(mObjectiveFunctionCalls == 0)
    force_update = true;
  
  mForceUpdate = mForceUpdate || force_update;
  mObjectiveFunctionCalls++;

// TODO: remove
  mForceUpdate = true;

  DoubleVector LLVector;
  if(by_row) LLVector = DoubleVector(mDatapoints.size());
  else       LLVector.push_back(0);

  // initialize coefficients
  CCoefs coefs =  CCoefs(raw_coefs, CCoefs::FnConstrain, false, mCoefConstraints);

  for(size_t i=0; i < mUniqueDMRows.size(); i++)
  {
    if(mEnabled[i] ) 
    {
        _dbg((0, "updated parameter for dm row %d/%d", i, mUniqueDMRows.size() )); 
        CDesignMatrixRow& cur_dmrow = mUniqueDMRows[i];
        std::vector<int>& datapoint_idxs = cur_dmrow.DatapointIndices();
        bool params_changed = cur_dmrow.UpdateParameters(coefs, mForceUpdate);
        _dbg((0, "updated parameters"));
        for(size_t j=0; j < datapoint_idxs.size(); j++)
        {
          int idx = datapoint_idxs[j];
        _dbg((0, "updating datapoint %d", idx));
          CDataPoint* cur_datapoint = mDatapoints[idx];
        _dbg((0, "updating datapoint %d (0x%x)", idx, cur_datapoint));
          if(params_changed)
            cur_datapoint->UpdateParameters(cur_dmrow, tolerate_imprecisions);
          double cur_LL = cur_datapoint->LogLik();
          save_LL(LLVector, cur_LL, idx, by_row);
        }
    }
  }

  mLastCoefs = coefs;

/*
  // this is mainly for debugging purposes
  // reset dprime and criterion
  for(size_t i=0; i < mDatapoints.size(); i++) {
    if( mEnabled[i] ) {
        mDatapoints[i].ResetDprime();
        mDatapoints[i].ResetCriterion();
    }
  }
*/  
  return LLVector;
}

Rcpp::DoubleVector CDataContainer::ObjectiveFunctionGradient(DoubleVector& raw_coefs, bool by_row, bool tolerate_imprecisions)
{
  _dbg_method("ObjectiveFunctionGradient", 1);

  DoubleVector gradient(raw_coefs.size(), 0.0);
  
/*
// initialize coefficients
  CCoefs coefs =  CCoefs(raw_coefs, CCoefs::FnConstrain, false, mCoefConstraints);

  if( !mDM.HasParameter( parameter_corr_mrsat ) )  {
    for(size_t i=0; i < mDatapoints.size(); i++) {
      if( mEnabled[i] ) {
          gradient += mDatapoints[i].ComputeLogLikGradient(coefs);
      }
    }

  } else {
      DoubleVector invalid_gradient(raw_coefs.size(), R_NaN);
      return invalid_gradient;
  }
  return gradient;
  */
      DoubleVector invalid_gradient(raw_coefs.size(), R_NaN);
      return invalid_gradient;
}


void CDataContainer::DetermineZeroRows(LogicalVector& zero_columns, std::vector<bool>& row_selected, bool all)
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
      int col_idx = mCoefConstraints.FindCoefIndex(column_name);

      for(size_t row_idx=0; row_idx < mUniqueDMRows.size(); row_idx++) 
      {
        double cur_value = mUniqueDMRows[row_idx][col_idx];
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
   
void CDataContainer::SelectSubset(LogicalVector& columns_zero, bool all) {
  DetermineZeroRows(columns_zero, mEnabled, all);
  mForceUpdate = true;
}

void CDataContainer::ResetSubset()
{
  for(size_t i=0; i < mEnabled.size(); i++)
    mEnabled[i] = true;
  mForceUpdate = true;
}

std::vector<int> CDataContainer::ReturnSelection()
{
  std::vector<int> selection;  
  for(size_t i=0; i < mUniqueDMRows.size(); i++) {
    if(mEnabled[i] ) {
      CDesignMatrixRow& cur_dmrow = mUniqueDMRows[i];
      std::vector<int>& datapoint_idxs = cur_dmrow.DatapointIndices();
      selection.insert(selection.end(), datapoint_idxs.begin(), datapoint_idxs.end());
    }
  }
  return selection;
}

