#ifndef HMM_H_
#define HMM_H_
#include <math.h>
#include <boost/function.hpp>
#include <Eigen/Dense>
#include "analysis/AnalysisBase.h"
#include "common/Track.h"

class HMM : public Analysis<float, int>
{
public: 
  typedef Eigen::ArrayXf VectorType;
  typedef Eigen::ArrayXXf MatrixType;
  typedef Eigen::ArrayXi StateVectorType;
  typedef Eigen::ArrayXXi StateMatrixType;
  
  using Analysis<float, int>::TrackInPtr;
  using Analysis<float, int>::TrackOutPtr;
  
  HMM() 
  : rand_init_probs_(true)
  , rand_trans_probs_(true)
  , has_trans_prior_(false)
  , has_init_prior_(false)
  {}
  
  virtual ~HMM() {}
 
  void set_transition(const MatrixType& m)
  { trans_ = m; rand_trans_probs_ = false; }

  const MatrixType& transition() const
  { return trans_; }

  void set_init_probs(const VectorType& init)
  { init_ = init; rand_init_probs_ = false; }

  const VectorType& init_probs() const 
  { return init_; }

  void set_init_probs_prior(const VectorType& prior)
  { init_prior_ = prior; has_init_prior_ = true; }

  const VectorType& init_probs_prior() const
  { return init_prior_; }

  // pseudo counts for transition matrix
  void set_trans_prior(const MatrixType& prior)
  { trans_prior_ = prior; has_trans_prior_ = true; }

  const MatrixType& trans_prior() const
  { return trans_prior_; }

  int num_states() const 
  { return num_states_; }

  void set_num_states(int s) 
  { num_states_ = s; }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  void ViterbiDecode(StateVectorType& output);
  float LogProb();
  
protected:
  void FitEM();
  void FwdBack();

  MatrixType CountTransitions(const StateVectorType& curr_states);

  // Compute fits the model
  virtual void ComputeAnalysis(TrackOutPtr output) = 0;
  
  // Implemented in the derived class b/c depends on the 
  // Particular emission distribution for the model.
  virtual void UpdateSoftEvidence() = 0;
  
  using Analysis<float, int>::input_;
  
  bool rand_init_probs_;
  bool rand_trans_probs_;
  bool has_trans_prior_;
  bool has_init_prior_;
  
  StateVectorType curr_states_;
  MatrixType soft_evidence_;
  
  MatrixType trans_; 
  MatrixType trans_prior_;
  VectorType init_;
  VectorType init_prior_;
  
  int num_states_;
};


class GaussHMM : public HMM
{
public:  
  using HMM::TrackInPtr;
  using HMM::TrackOutPtr;

  GaussHMM() {}
  virtual ~GaussHMM() {}
  
  virtual void ComputeAnalysis(TrackOutPtr output);
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
  
};




#endif
