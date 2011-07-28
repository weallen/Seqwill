#include "analysis/HMM.h"

float 
HMM::LogProb(const MatrixType& softev)
{
  float logprob = 0.0;
  for (int i = 0; i < softev.cols(); ++i) {
    logprob += log(softev.col(i).maxCoeff());
  }
  return logprob;
}

void
HMM::ViterbiDecode(const VectorType& init, const MatrixType& trans, const MatrixType& softev, StateVectorType& path)
{ 
  int T = softev.cols() - 1;
  path.resize(T);
  // Delta holds the normalized values
  MatrixType delta = MatrixType::Zero(num_states_, T);
  
  // psi stores the indices of the maximum value for each
  StateMatrixType psi = StateMatrixType::Zero(num_states_, T);
  
  // Initialize last element as soft evidence weighted
  // By prior probabilities
  delta.col(T) = init * softev.col(T);  
  VectorType v;
  MatrixType::Index idx;
  for (int t=1; t <= T; ++t) {
    for (int j=0; j < num_states_; ++j) {
      // find the unnormalized values
      v = delta.col(t-1) * trans.col(j);
      delta(j,t) = v.maxCoeff(&idx) * softev(j,t);
      psi(j,t) = (int) idx;
    }
    // Normalize 
    delta.col(t) /= delta.col(t).sum();
  }
  
  // Traceback
  // First find the 

  delta.col(T).maxCoeff(&idx);
  path(T) = (int)idx;  
  for (int t = T - 1; t >= 0; --t) {
    path(t) = psi(path(t+1), t+1);
  }
}


HMM::MatrixType 
HMM::CountTransitions(const StateVectorType& curr_states)
{
  MatrixType m = MatrixType::Zero(num_states_, num_states_);
  
  // Count transitions
  int curr_state;
  int next_state;
  for (int pos = 0; pos < curr_states.size()-1; ++pos) {
    curr_state = curr_states(pos);
    next_state = curr_states(pos+1);
    m(curr_state, next_state) += 1;
  }
  // Normalize matrix
  for (int i = 0; i < num_states_; ++i) {
    for (int j = 0; j < num_states_; ++j) {
      m(i,j) = m(i,j) / m.row(i).sum();
    }
  }
  return m;
}

HMM::MatrixType 
HMM::TwoSliceSum(const MatrixType& transmat, const MatrixType& softev,
                 const MatrixType& alpha, const MatrixType& beta)
{
  MatrixType xi_summed = MatrixType::Zero(num_states_, num_states_);
  int T = softev.cols();
  VectorType b;
  MatrixType xit;
  
  for (int t = T-2; t > 0; --t) {
    b = beta.col(t+1) * softev.col(t+1);
    xit = transmat * (alpha.col(t).matrix() * b.matrix().transpose()).array();
    // XXX This might not work: equiv to xi_summed = xi_summed + xit ./ sum(xit(:))?
    xi_summed += xit / xit.colwise().sum();
  }
  return xi_summed;
}

float 
HMM::FwdBack(const MatrixType& transmat, const VectorType& init, 
             const MatrixType& softev, MatrixType& alpha, MatrixType& beta, 
             MatrixType& gamma)
{  
  int T = softev.cols();
  if (gamma.cols() != T && gamma.rows() != num_states_) {
    gamma.resize(num_states_, T);
  }
  float loglik;
  
  DEBUGLOG("Filter forward");
  FilterFwd(transmat, softev, init, loglik, alpha);
  
  DEBUGLOG("Smooth backward");
  SmoothBack(transmat, softev, beta);
  
  gamma = alpha * beta;
  for (int t = 0; t < T; ++t) {
    gamma.col(t) /= gamma.col(t).sum();
  }
  return loglik;
}

void 
HMM::FilterFwd(const MatrixType& transmat, const MatrixType& softev, 
               const VectorType& init, float& loglik, MatrixType& alpha)
{
  int T = (int) softev.cols();
  
  if (alpha.cols() != T && alpha.rows() != num_states_) {
    alpha.resize(num_states_, T);
  }
  VectorType scale = VectorType::Zero(T);

  alpha.col(1) = init * softev.col(1);
  scale(1) = alpha.col(1).sum();
  alpha.col(1) /= scale(1);
  Eigen::MatrixXf at = transmat.matrix().transpose();
  for (int t = 1; t < T; ++t) {
    alpha.col(t) = (transmat.matrix() * alpha.col(t-1).matrix()).array();
    alpha.col(t) *= softev.col(t);
    scale(t) = alpha.col(t).sum();
    alpha.col(t) /= scale(t);
  }
  loglik = scale.log().sum();
}

void 
HMM::SmoothBack(const MatrixType& transmat, const MatrixType& softev, 
                MatrixType& beta)
{
  int K = (int)softev.rows();
  int T = (int)softev.cols();
  
  beta.resize(K, T);  
  for (int t = T-2; t > 0; --t) {    
    // beta(:,t) = trans_ * (beta(:,t+1) .* soft_evidence_(:,t+1))
    beta.col(t) = transmat.matrix() * (beta.col(t+1) * softev.col(t+1)).matrix();
    // normalize
    beta.col(t) /= beta.col(t).sum();
  }
}

void 
HMM::FitEM() 
{
  
  float loglik;
  int length = input_->size();
  MatrixType bi;
  MatrixType xi_summed;
  MatrixType alpha(num_states_, num_states_);
  MatrixType beta(num_states_, 2);
  MatrixType gamma;
  MatrixType trans = trans_prior_;
  VectorType init = init_prior_;
  VectorType start_counts = VectorType::Zero(num_states_, num_states_);

  for (int n=0; n < 100; ++n) {
    DEBUGLOG("EM step " + Stringify(n));
    //
    // E STEP  
    //
    for (int i=0; i < length; ++i) {
      // XXX Not sure if this is correct.
      bi = soft_evidence_.block(num_states_, 2,1, i);
      loglik += FwdBack(trans, init, bi, alpha, beta, gamma);
      xi_summed = TwoSliceSum(trans, bi, alpha, beta);
    }
    
    //
    // M STEP
    // 
    
    UpdateEmissionDist();
  }
}


//----------------------------------------------------------------

void
GaussHMM::ComputeAnalysis()
{
  DEBUGLOG("Fitting model by EM");
  FitEM();
  
  DEBUGLOG("Predicting most likely state sequence");
  StateVectorType v;
  ViterbiDecode(v);
}

void
GaussHMM::UpdateSoftEvidence()
{  
  soft_evidence_.resize(num_states_, input_->size());
  for (int i = 0; i < soft_evidence_.cols(); ++i) {
    for (int j = 0; j < soft_evidence_.rows(); ++j) {
      soft_evidence_(j, i) = emit_[j].pdf(input_->get(i));
    }
  }
}

void
GaussHMM::UpdateEmissionDist(const StateVectorType& curr_states)
{
  std::vector<float> obs;
  for (int i =0; i < num_states_; ++i) {
    obs.clear();
    for (int t = 0; t < curr_states.cols(); ++t) {
      if (curr_states(t) == i) {
        obs.push_back(curr_states(t));
      }
    }
    float mean = 0.0;
    float stddev = 0.0;
    for (size_t k = 0; k < obs.size(); ++k) {
      mean += obs[k];
    }
    mean /= (float) obs.size();
    for (size_t k = 0; k < obs.size(); ++k) {
      stddev += (obs[k] - mean) * (obs[k] - mean);
    }
    stddev /= (float) obs.size();
    emit_[i].set_mean(mean);
    emit_[i].set_stddev(stddev);
  }    
}