#ifndef HMM_H_
#define HMM_H_
#include <math.h>
#include <boost/function.hpp>
#include <Eigen/Dense>
#include "analysis/AnalysisBase.h"
#include "analysis/Dist.h"
#include "common/Track.h"
#include "base/Log.h"
#include "base/StringUtil.h"
#include "analysis/Random.h"

class HMM : public Analysis<float, int>
{
public: 
  typedef Eigen::ArrayXd VectorType;
  typedef Eigen::ArrayXXd MatrixType;
  typedef Eigen::ArrayXi StateVectorType;
  typedef Eigen::ArrayXXi StateMatrixType;
  
  using Analysis<float, int>::TrackInPtr;
  using Analysis<float, int>::TrackOutPtr;

  HMM() 
  : rand_init_probs_(true)
  , rand_trans_probs_(true)
  , has_trans_prior_(false)  
  , has_init_prior_(false)
  , num_states_(1)
  {  }
  
  HMM(int num_states)
  : rand_init_probs_(true)
  , rand_trans_probs_(true)
  , has_trans_prior_(false)  
  , has_init_prior_(false)
  , num_states_(num_states)
  {  }

  virtual ~HMM() {}
 
  void Init();
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
  { num_states_ = s; NumStatesChanged(); }

  // OUTPUT: log(p(y(1:T))
  float LogProb(const MatrixType& softev);  
  
  void FitEM();
  void FitBlockedGibbs();
  
  void Decode(StateVectorType& path);
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  
protected:  
  virtual void NumStatesChanged() = 0;
  
  // OUTPUT: 
  void ViterbiDecode(const VectorType& init, const MatrixType& trans,
                     const MatrixType& softev, StateVectorType& path);
  
  
  // Computes p(S(t) = i | y(1:T))
  // OUTPUT: gamma(i,t) = p(S(t) = i | y(1:T)), returns loglik
  // The probability of each hidden variable being in a particular 
  // state given all the evidence.
  // gamma = alpha * beta
  // i.e. p(S(t) = i | y(1:T)) = p(S(t) = i | y(1:t), y(t+1:T))
  //      \propto p(S(t) = i | y(1:t)) * p(y(t+1:T) | S(t) = i)  
  float FwdBack(const MatrixType& transmat, const VectorType& init, 
                const MatrixType& softev, MatrixType& alpha, MatrixType& beta, 
                MatrixType& gamma);

  // Computes p(S(t) = i | y(1:t))
  void FilterFwd(const MatrixType& transmat, const MatrixType& softev, 
                 const VectorType& init, float& loglik, MatrixType& alpha);
  
  void SampleBack(const MatrixType& transmat, const MatrixType& softev, 
                  const MatrixType& alpha, StateVectorType& curr_states);
  
  // Computes p(y(t+1:T) | S(t) = i)
  void SmoothBack(const MatrixType& transmat, const MatrixType& softev, 
                  MatrixType& beta);
  
  // OUTPUT: new trans_ matrix
  MatrixType CountTransitions(const StateVectorType& curr_states);

  // Computes the sum of two-slice distribution over hidden states
  // INPUT: alpha and beta from FwdBack
  // Also uses current soft_evidence_ and trans_
  // OUTPUT: xi_summed(i,j) = \sum_t=2^T p(S(t) = i, S(t+1) = j | y(1:T))
  // Expected sufficient statistics for the transition matrix, 
  // for a given observation sequence.
  MatrixType TwoSliceSum(const MatrixType& transmat, const MatrixType& softev,
                         const MatrixType& alpha, const MatrixType& beta);
  
  // Compute fits the model
  virtual void ComputeAnalysis() = 0;
  
  // Implemented in the derived class b/c depends on the 
  // Particular emission distribution for the model.
  // OUTPUT: softev(i,t) = p(y(t) | S(t) = i)
  virtual void UpdateSoftEvidence(MatrixType& softev) = 0;
  virtual void UpdateEmissionDistEM(const MatrixType& weights) = 0;
  virtual void UpdateEmissionDistGibbs(Rng& r, const StateVectorType& states) = 0;

  using Analysis<float, int>::input_;
  
  bool rand_init_probs_;
  bool rand_trans_probs_;
  bool has_trans_prior_;
  bool has_init_prior_;
  
  // soft_evidence_(i,t) = p(y(t) | S(t) = i)
  // Moved to inside FitEM
  //MatrixType soft_evidence_;
  
  // trans_(i,j) = p(S(t) = j | S(t-1) = i)
  MatrixType trans_; 
  
  // trans_prior_(i,j) = p(S(t) = j, S(t-1) = i | alpha)
  MatrixType trans_prior_;
  
  // init_(i) = p(S(i) = i)
  VectorType init_;
  
  // init_prior_(i) = p(S(i) = i | alpha)
  VectorType init_prior_;
  
  int num_states_;
  
  friend class HMMTest;
};

class GaussHMM : public HMM
{
public:  
  using HMM::TrackInPtr;
  using HMM::TrackOutPtr;
  using Analysis<float, int>::Compute;
  
  GaussHMM();
  GaussHMM(int num_states);
  virtual ~GaussHMM() {}
  
  virtual void ComputeAnalysis();
  virtual void UpdateSoftEvidence(MatrixType& softev);
  virtual void Compute() { ComputeAnalysis(); }
  virtual void UpdateEmissionDistEM(const MatrixType& weights);
  virtual void UpdateEmissionDistGibbs(Rng& r, const StateVectorType& states);
  virtual void NumStatesChanged();
    
  void set_emit(const std::vector<GaussDist>& emits)
  { emit_ = emits; }
  
  const std::vector<GaussDist>& emit() const { return emit_; }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
protected:
  ChiSqDist var_prior;
  StudentTDist mean_prior;
  std::vector<GaussDist> emit_;
};


class GaussMultiTrackHMM : public GaussHMM
{
public:
  using HMM::TrackInPtr;
  using HMM::TrackOutPtr;
  using Analysis<float, int>::Compute;
  
  GaussMultiTrackHMM();
  GaussMultiTrackHMM(int num_states);
  virtual ~GaussMultiTrackHMM() {}
  
  virtual void UpdateSoftEvidence(MatrixType& softev);
  virtual void UpdateEmissionDistEM(const MatrixType& weights);
  virtual void UpdateEmissionDistGibbs(Rng& r, const StateVectorType& states);
  virtual void NumStatesChanged();
  
  void set_emit(const std::vector<MVGaussDist>& emits)
  { emits_ = emits; }    
  const std::vector<MVGaussDist>& emit() const { return emits_; }
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  virtual void set_input(TrackInPtr input)
  { input_ = input; tracks_[0] = input; }
  
  virtual void add_track(TrackInPtr input);
  
protected:
  
  std::vector<TrackInPtr> tracks_;
  std::vector<MVGaussDist > emits_;
};






#endif
