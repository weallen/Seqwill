#include "analysis/HMM.h"

void
HMM::Init()
{
  if (rand_init_probs_) {
    DEBUGLOG("Initializing HMM with random pi");    
    init_ = VectorType::Random(num_states_).abs();    
    init_ /= init_.sum();
  } if (rand_trans_probs_) {
    DEBUGLOG("Initializing HMM with random trans probs");
    trans_ = MatrixType::Random(num_states_, num_states_).abs();
    for (int i = 0; i < trans_.cols(); ++i) {
      trans_.col(i) /= trans_.col(i).sum();
    }
  }
  if (!has_trans_prior_) {
    DEBUGLOG("Initializing HMM with no transition prior");
    trans_prior_ = MatrixType::Zero(num_states_, num_states_);
  }
  if (!has_init_prior_) {
    DEBUGLOG("Initializing HMM with no pi prior");
    init_prior_ = MatrixType::Zero(num_states_, num_states_);
  }
}

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
  int T = softev.cols();
  path.resize(T);
  // Delta holds the normalized values
  MatrixType delta = MatrixType::Zero(num_states_, T);
  
  // psi stores the indices of the maximum value for each
  StateMatrixType psi = StateMatrixType::Zero(num_states_, T);
  
  // Initialize last element as soft evidence weighted
  // By prior probabilities
  delta.col(T-1) = init * softev.col(T-1);  
  VectorType v;
  MatrixType::Index idx;
  for (int t=1; t < T; ++t) {
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

  delta.col(T-1).maxCoeff(&idx);
  path(T-1) = (int)idx;  
  for (int t = T - 2; t >= 0; --t) {
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
  /*
  for (int i = 0; i < num_states_; ++i) {
    for (int j = 0; j < num_states_; ++j) {
      m(i,j) = m(i,j) / m.row(i).sum();
    }
  }*/
  return m;
}

// Outputs the sufficient statistics for the E step
HMM::MatrixType 
HMM::TwoSliceSum(const MatrixType& transmat, const MatrixType& softev,
                 const MatrixType& alpha, const MatrixType& beta)
{
  MatrixType xi_summed = MatrixType::Zero(num_states_, num_states_);
  int T = softev.cols();
  VectorType b;
  MatrixType xit;
  
  for (int t = T-2; t >= 0; --t) {
    // multiply by soft evidence
    b = beta.col(t+1) * softev.col(t+1);
    
    xit = transmat * (alpha.col(t).matrix() * b.matrix().transpose()).array();    
    xi_summed += xit / xit.sum();
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
  

  //DEBUGLOG("Filter forward");
  FilterFwd(transmat, softev, init, loglik, alpha);
  

  //DEBUGLOG("Smooth backward");
  SmoothBack(transmat, softev, beta);
  int i = 0; int j = 0;
  beta.minCoeff(&i, &j);
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
  alpha.col(0) = init * softev.col(0);
  scale(0) = alpha.col(0).sum();
  alpha.col(0) /= scale(0);
  Eigen::MatrixXf at = transmat.matrix().transpose();
  for (int t = 1; t < T; ++t) {
    alpha.col(t) = (at.matrix() * alpha.col(t-1).matrix()).array();
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
  if (beta.cols() != T && beta.rows() != K) 
    beta.resize(K, T);  
  for (int k = 1; k < num_states_; ++k) 
    beta(k, T-1) = 1.0;
  for (int t = T-2; t >= 0; --t) {    
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
  MatrixType alpha;
  MatrixType beta;
  MatrixType gamma;
  MatrixType trans = trans_;
  VectorType init = init_;
  VectorType start_counts = VectorType::Zero(num_states_);
  MatrixType trans_counts = MatrixType::Zero(num_states_, num_states_);
  VectorType weights = VectorType::Zero(num_states_);
  MatrixType softev = MatrixType::Zero(num_states_, length);
  MatrixType bi;
  MatrixType xi_summed;
  
  for (int n=0; n < 20; ++n) {
    DEBUGLOG("EM step " + Stringify(n));
    //
    // E STEP  
    //
    UpdateSoftEvidence(softev);

    
    loglik = FwdBack(trans, init, softev, alpha, beta, gamma);
    trans_counts = TwoSliceSum(trans, softev, alpha, beta);
    start_counts = gamma.col(0);
    std::cerr << "LOGLIK: " << loglik << std::endl;
    
    //
    // M STEP
    // 
    //ViterbiDecode(init, trans, gamma, curr_path);
    //trans = CountTransitions(curr_path) + trans_prior_;
    //for (int i = 0; i < trans.cols(); ++i) {
    //  trans.col(i) /= trans.col(i).sum();
    //}
    //for (int i = 0; i < num_states_; ++i) {
    //  init(i) = (curr_path == i).count() + init_prior_(i);
    //  DEBUGLOG("Predicts " + Stringify((curr_path == i).count()) + " in state " + Stringify(i));
    //}    
    //init /= init.sum();  
    
    init = start_counts + init_prior_;
    init /= init.sum();
    trans = trans_counts + trans_prior_;
    for (int i = 0; i < trans.cols(); ++i) {
      trans.col(i) /= trans.col(i).sum();
    }
    UpdateEmissionDist(gamma);
  }
  trans_ = trans;
  init_ = init;  
}


//----------------------------------------------------------------


GaussHMM::GaussHMM()
: HMM()
, emit_(1)
{}

GaussHMM::GaussHMM(int num_states)
: HMM(num_states)
, emit_(num_states)
{}

void
GaussHMM::NumStatesChanged()
{
  emit_.resize(num_states_);
}

void
GaussHMM::ComputeAnalysis()
{
  DEBUGLOG("Fitting model by EM");
  FitEM();
  
  DEBUGLOG("Predicting most likely state sequence");
  StateVectorType v;
  //ViterbiDecode(v);
}

void
GaussHMM::UpdateSoftEvidence(MatrixType& softev)
{  
  if (softev.cols() != input_->size() || softev.rows() != num_states_) 
    softev.resize(num_states_, input_->size());
  for (int i = 0; i < softev.cols(); ++i) {
    for (int j = 0; j < softev.rows(); ++j) {
      softev(j, i) = emit_[j].pdf(input_->get(i));
    }
  }
}

void
GaussHMM::UpdateEmissionDist(const MatrixType& weights)
{
  for (int k = 0; k < num_states_; ++k) {
    float norm = weights.row(k).sum();
    float mean = 0.0;
    float stddev = 0.0;      
    for (size_t i = 0; i < input_->size(); ++i) {
      mean += input_->get(i);
    }
    mean /= norm;
    for (size_t i = 0; i < input_->size(); ++i) {
      stddev += weights(k,i) * (input_->get(i) - mean) * (input_->get(i) - mean);
    }
    stddev /= norm;
    emit_[k].set_mean(mean);
    emit_[k].set_stddev(sqrt(stddev));
  }      
}