#ifndef HMM_H_
#define HMM_H_
#include <boost/function.hpp>
#include <Eigen/Dense>

#include "analysis/AnalysisBase.h"

template <class EmitDistT>
class HMM : public Analysis<float, float>
{
public: 
  typedef Eigen::Vector<float, Eigen::Dynamic, 1> VectorType;
  typedef Eigen::Vector<float Nstates, Nstates> MatrixType;

  HMM() 
  : rand_init_probs_(true)
  , rand_trans_probs_(true)
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
  { init_prior_ = prior; }

  const VectorType& init_probs_prior() const
  { return init_prior_; }

  // pseudo counts for transition matrix
  void set_trans_prior(const MatrixType& prior)
  { trans_prior_ = prior; }

  const MatrixType& trans_prior() const
  { return trans_prior_; }

  void set_emit_fn(const EmitFnType& fun) 
  { emit_fn_ = fun; }

  const EmitFnType& emit_fn() const
  { return emit_fn_; }

  int num_states() const 
  { return num_states_; }

  void set_num_states(int s) 
  { num_states_ = s; }

  float LogProb();
  void ViterbiDecode();
  
private:
  void FitEM();
  void FwdBack(const VectorType& pi, 
               const MatrixType& A,
               Eigen::MatrixXf* gamma,
               Eigne::MatrixXf* alpha,
               float* logp);
   
  void CountTransitions(MatrixType* out);

  // Compute fits the model
  virtual void Compute();

  bool rand_init_probs_;
  bool rand_trans_probs_;
  MatrixType trans_; 
  MatrixType trans_prior_;
  VectorType init_;
  VectorType init_prior_;
  int num_states_;
  std::vector<EmitDistT> emit_dists_;
};

template<>
void
HMM<GaussDist>::UpdateEmitEM()
{
  
}

template <class EmitTypeT>
void 
HMM<EmitTypeT>::FitEM() 
{
  MatrixType start_counts;
  MatrixType trans_counts;
  int len = input_.size();
  float logp;
  MatrixXf alpha(num_states(), len); 
  MatrixXf beta(num_states(), len);
  for (int i = 0; i < len; ++i) {
    FwdBack( &logp);
    DEBUGLOG(logp);
  }
}

template <class EmitTypeT>
void 
HMM<EmitTypeT>::FwdBack()
{
}

template <class EmitTypeT>
void 
HMM<EmitTypeT>::ViterbiDecode(std::vector<int>& output)
{
}

template <class EmitTypeT>
MatrixType 
HMM<EmitTypeT>::CountTransitions()
{
  return m;
}



#endif
