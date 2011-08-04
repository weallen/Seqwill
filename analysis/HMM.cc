#include "analysis/HMM.h"

void
HMM::Init()
{
  if (rand_init_probs_) {
    init_ = VectorType::Random(num_states_).abs();    
    init_ /= init_.sum();
        DEBUGLOG("Initializing HMM with random pi");    
  } if (rand_trans_probs_) {
    DEBUGLOG("Initializing HMM with random trans probs");
    trans_ = MatrixType::Random(num_states_, num_states_).abs();
    for (int i = 0; i < trans_.rows(); ++i) {
      trans_.row(i) /= trans_.row(i).sum();
    }
  }
  if (!has_trans_prior_) {
    DEBUGLOG("Initializing HMM with no transition prior");
    trans_prior_ = MatrixType::Zero(num_states_, num_states_);
  }
  if (!has_init_prior_) {
    DEBUGLOG("Initializing HMM with no pi prior");
    init_prior_ = VectorType::Zero(num_states_);
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
HMM::Decode(StateVectorType& path)
{ 
  MatrixType softev; 
  UpdateSoftEvidence(softev);
  ViterbiDecode(init_, trans_, softev, path); 
}

void
HMM::ViterbiDecode(const VectorType& init, const MatrixType& trans, 
                   const MatrixType& softev, StateVectorType& path)
{ 
  int T = softev.cols();
  if (path.size() != T) {
    path.resize(T);
  }
  // Delta holds the normalized values
  MatrixType delta = MatrixType::Zero(num_states_, T);
  
  // psi stores the indices of the maximum value for each
  StateMatrixType psi = StateMatrixType::Zero(num_states_, T);
  
  // Initialize last element as soft evidence weighted
  // By prior probabilities
  delta.col(0) = init * softev.col(0);  
  delta.col(0) /= delta.col(0).sum();
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

  for (int pos = 1; pos < curr_states.size(); ++pos) {
    m(curr_states(pos-1), curr_states(pos)) += 1;
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

  // Normalize
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
  Eigen::MatrixXd at = transmat.matrix().transpose();
  
  alpha.col(0) = init * softev.col(0);
  scale(0) = alpha.col(0).sum();
  alpha.col(0) /= scale(0);  

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
HMM::FitBlockedGibbs()
{
  int length = input_->size();
  Rng r;
  MultiDist m;
  
  // Randomly initialize states
  m.set_vals(VectorType::Constant(num_states_, 1.0/((float) num_states_)));
  StateVectorType curr_states = StateVectorType::Zero(length);
  for (int i = 0; i < length; ++i) {
    curr_states(i) = m.Sample(r.rng());
  }
  MatrixType alpha;
  MatrixType transmat = trans_;
  MatrixType trans_counts;
  MatrixType init = init_;
  MatrixType softev;
  DirichletDist trans_dist;
  
  float loglik;
  for (int n = 0; n < 100; ++n) {
    DEBUGLOG("Blocked Gibbs step: " + Stringify(n));
    UpdateSoftEvidence(softev);
    std::cerr << "Filter fwd" << std::endl;
    FilterFwd(transmat, softev, init, loglik, alpha);
    
    
    std::cerr << "Sample back" << std::endl;
    SampleBack(transmat, softev, alpha, curr_states);
    
    // Sample transition matrix, regularized by adding a pseudocount 
    // XXX Make pseudocount a parameter of the model
    trans_counts = CountTransitions(curr_states);
    for (int i = 0; i < num_states_; ++i) {
      trans_dist.set_alpha(trans_counts.row(i) + 1);
      transmat.row(i) = trans_dist.Sample(r.rng());
    }
    std::cerr << "TRANSMAT: " << transmat << std::endl;
    std::cerr << (curr_states == 0).count() << std::endl;
    std::cerr << (curr_states == 1).count() << std::endl;
    // Sample emission parameters
    std::cerr << "Update emission" << std::endl;
    UpdateEmissionDistGibbs(r, curr_states);
  }
  trans_ = transmat;
  init_ = init;
}

void
HMM::SampleBack(const MatrixType& transmat, const MatrixType& softev, 
                const MatrixType& alpha, StateVectorType& curr_states)
{
  MultiDist multi(num_states_);
  Rng r;
  VectorType gamma;

  int K = (int)softev.rows();
  int T = (int)softev.cols();
  MatrixType beta(K,T);
  
  for (int k = 1; k < num_states_; ++k) 
    beta(k, T-1) = 1.0;
  gamma = beta.col(T-1) * alpha.col(T-1);
  multi.set_vals(gamma);
  curr_states(T-1) = multi.Sample(r.rng());
  for (int t = T-2; t >= 0; --t) {    
    // beta(:,t) = trans_ * (beta(:,t+1) .* soft_evidence_(:,t+1))
    beta.col(t) = transmat.matrix() * (beta.col(t+1) * softev.col(t+1)).matrix();    
    // normalize
    beta.col(t) /= beta.col(t).sum();
    gamma = beta.col(t) * alpha.col(t);
    multi.set_vals(gamma);
    curr_states(t) = multi.Sample(r.rng());
  }

/*  int T = alpha.cols();
  
  if (curr_states.size() != alpha.cols()) {
    curr_states.resize(alpha.cols());
  }
  Eigen::VectorXd curr_vals = Eigen::VectorXd::Zero(num_states_);
  // sample z_T ~ p(z_T = i| x_{1:T}) = \alpha_T(i)
  multi.set_vals(alpha.col(T-1));
  curr_states(T-1) = multi.Sample(r.rng());
  // sample z_t^* ~ p(z_t|z_{t+1}^*, x_{1:T}) 
  for (int t = T-2; t >= 0; --t) {
    int j = curr_states(t+1);
    for (int i = 0; i < num_states_; ++i) {      
      curr_vals(i) = (trans_(i,j) * alpha(i,t) * softev(j, t+1)) / alpha(j,t+1);
    }
    multi.set_vals(curr_vals);
    curr_states(t) = multi.Sample(r.rng());
  }  */
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

  StateVectorType v;
  float delta_loglik = INFINITY;
  float old_loglik = 0.0;
  int n = 0;
  while (delta_loglik > 0.001) {
    n++;
    DEBUGLOG("EM step " + Stringify(n));
    //
    // E STEP  
    //
    UpdateSoftEvidence(softev);
    
    loglik = FwdBack(trans, init, softev, alpha, beta, gamma);    
    trans_counts = TwoSliceSum(trans, softev, alpha, beta);

    std::cerr << "LOGLIK: " << loglik << std::endl;
    //
    // M STEP
    // 
    
    init = gamma.col(1) + init_prior_;
    init /= init.sum();
    trans = trans_counts + trans_prior_;
    VectorType z = trans.rowwise().sum();
    for (int i = 0; i < num_states_; ++i) {
      trans.row(i) /= z; 
    }
    UpdateEmissionDistEM(gamma);
    delta_loglik = abs(old_loglik - loglik);
    old_loglik = loglik;
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
}

void
GaussHMM::UpdateSoftEvidence(MatrixType& softev)
{  
  if (softev.cols() != (int)input_->size() || softev.rows() != num_states_) 
    softev.resize(num_states_, input_->size());
  for (int i = 0; i < softev.cols(); ++i) {
    for (int j = 0; j < softev.rows(); ++j) {
      
      // XXX HACK ALERT -- ADD MACHINE EPSILON TO PREVENT UNDERFLOW
      softev(j, i) = emit_[j].Pdf((double)input_->get(i)) + 1.11e-16;
    }
  }
}

void
GaussHMM::UpdateEmissionDistEM(const MatrixType& weights)
{
  for (int k = 0; k < num_states_; ++k) {
    double norm = weights.row(k).sum();
    double mean = 0.0;
    double stddev = 0.0;      
    for (size_t i = 0; i < input_->size(); ++i) {
      mean += weights(k,i) * (double)input_->get(i);
    }
    mean /= norm;
    for (size_t i = 0; i < input_->size(); ++i) {
      stddev += weights(k,i) * ((double)input_->get(i) - mean) * ((double)input_->get(i) - mean);
    }
    stddev /= norm;    
    emit_[k].set_mean(mean);
    emit_[k].set_stddev(sqrt(stddev));
    //    std::cerr << "State " << k << " mean " << mean << " std " << stddev << std::endl;
  }      
}


// From paper "Inferring a Gaussian distribution" by T. Minka, 2001
// And pg 83 of Bayesian Data Analysis ("Samping from the Joint Posterior
// Distribution" of a normal distribution)
void
GaussHMM::UpdateEmissionDistGibbs(Rng& r, const StateVectorType& states)
{
  GaussDist mean_dist;
  InvChiSqDist var_dist;
  
  // sample parameters of Gaussian dist from posterior distribution conditioned on current
  // assignment of the latent variables.
  //
  
  for (int i = 0; i < num_states_; ++i) {
    int state_count = 0;
    float mean = 0.0;
    float var = 0.0;
    double mean_sum = 0.0;
    double var_sum = 0.0;
    double tau = 0.0;
    
    // Hyperparams
    double s0 = 1.0;
    double m0 = 0.01;
    double t0 = 1.0;
    int nu = 1;
    
    for (size_t t = 0; t < input_->size(); ++t) {
      if (states(t) == i) {
        state_count++; 
        mean_sum += input_->get(t);
      }
    }
    
    mean = mean_sum / (double) state_count;
    for (int j = 0; j < 3; ++j) {
      // Sample each several times
      var_sum = 0.0;
      // Sample sigma_k | mu_k ~ InvChiSq(var, N)
      for (size_t t = 0; t < input_->size(); ++t) {
        if (states(t) == i) {
          var_sum += pow((input_->get(t) - mean), 2);
        }
      }
      var_dist.set_df(state_count - 1);
      var_dist.set_scale(var_sum / (state_count - 1));
      
      //  sample var from it's posterior
      var = var_dist.Sample(r.rng());
      
      // Sample mu_k | sigma_k ~ Gaussian(mean, var|data)
      mean_dist.set_stddev(1.0/(1.0/t0 + state_count / var));
      mean_dist.set_mean((m0 / t0 + 1.0 / var * mean_sum) / (1.0 / t0 + state_count/var));
      mean = mean_dist.Sample(r.rng());
    }
    emit_[i].set_stddev(var);  
    emit_[i].set_mean(mean);
    
    std::cerr << "State " << i << " mean " << emit_[i].mean() << " std " << emit_[i].stddev() << std::endl;
  }  
}



//-----------------------------------------------------------------

GaussMultiTrackHMM::GaussMultiTrackHMM()
: GaussHMM()
, emits_(1)
, tracks_(0)
{
}

GaussMultiTrackHMM::GaussMultiTrackHMM(int num_states)
: GaussHMM(num_states)
, emits_(num_states)
, tracks_(0)
{
}

void
GaussMultiTrackHMM::NumStatesChanged()
{
  emits_.resize(num_states_);
}

void
GaussMultiTrackHMM::add_track(GaussMultiTrackHMM::TrackInPtr input)
{
  if (tracks_.size() == 0) 
    input_ = input;
  tracks_.push_back(input);
  for (int i = 0; i < num_states_; ++i) {
    emits_[i] = std::vector<GaussDist>(tracks_.size()); 
  }
}

void 
GaussMultiTrackHMM::UpdateEmissionDistEM(const MatrixType& weights)
{
  int tnum = tracks_.size();
  Eigen::VectorXd v(tnum);
  std::vector<GaussDist> temp;  
  for (int k = 0; k < num_states_; ++k) {    
    double norm = weights.row(k).sum();
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(tnum);
    Eigen::VectorXd stddev = Eigen::VectorXd::Zero(tnum);
    for (size_t i = 0; i < input_->size(); ++i) {
      for (int tr = 0; tr < tracks_.size(); ++tr) {
        v(tr) = (double) tracks_[tr]->get(i);
      }
      mean += weights(k,i) * v;
    }
    mean /= norm;
    // make sure variance is a diagonal matrix -- assumes tracks are independent.
    for (size_t i = 0; i < input_->size(); ++i) {
      for (int tr = 0; tr < tracks_.size(); ++tr) {
	stddev(tr) += weights(k,i) * pow((double)tracks_[tr]->get(i) - mean(tr),2);
      }
    }
    stddev /= norm;    
    for (int tr = 0; tr < tracks_.size(); ++tr) {
      temp = emits_[k];
      temp[tr].set_mean(mean(tr));
      temp[tr].set_stddev(sqrt(stddev(tr)));
      std::cerr << "Track " << tr  << " Mean: " << temp[tr].mean() << " Stddev: " << temp[tr].stddev() << std::endl;
      emits_[k] = temp;
    }
  }        
}

void 
GaussMultiTrackHMM::UpdateEmissionDistGibbs(Rng& r, const StateVectorType& states)
{
  ERRORLOG("NOT IMPLMEENTED YET");
}

void
GaussMultiTrackHMM::UpdateSoftEvidence(MatrixType& softev)
{  
  if (softev.rows() != num_states_ || softev.cols() != (int)input_->size())
    softev.resize(num_states_, input_->size());
  int nt = tracks_.size();
  std::vector<GaussDist> temp;
  Eigen::VectorXd v(nt);
  for (int i = 0; i < softev.cols(); ++i) {
    for (int j = 0; j < softev.rows(); ++j) {
      softev(j, i) = 1.0;
      temp = emits_[j];
      for (int k = 0; k < nt; ++k) {
        softev(j,i) *= temp[k].Pdf((double)tracks_[k]->get(i)) + 1.11e-16;
      }
    }
  }
}






//-----------------------------------------------------------------

MVGaussMultiTrackHMM::MVGaussMultiTrackHMM()
: GaussHMM()
{
}

MVGaussMultiTrackHMM::MVGaussMultiTrackHMM(int num_states)
: GaussHMM(num_states)
{
}

void
MVGaussMultiTrackHMM::NumStatesChanged()
{
  emits_.resize(num_states_);
}

void
MVGaussMultiTrackHMM::add_track(GaussMultiTrackHMM::TrackInPtr input)
{
  if (tracks_.size() == 0) 
    input_ = input;
  tracks_.push_back(input);
}

void 
MVGaussMultiTrackHMM::UpdateEmissionDistEM(const MatrixType& weights)
{
  int tnum = tracks_.size();
  double temp;
  Eigen::VectorXd v(tnum);
  Eigen::MatrixXd ident = Eigen::MatrixXd::Identity(tnum, tnum);
  for (int k = 0; k < num_states_; ++k) {    
    double norm = weights.row(k).sum();
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(tnum);
    Eigen::MatrixXd stddev = Eigen::MatrixXd::Zero(tnum, tnum);
    for (size_t i = 0; i < input_->size(); ++i) {
      for (int tr = 0; tr < tracks_.size(); ++tr) {
        v(tr) = (double) tracks_[tr]->get(i);
      }
      mean += weights(k,i) * v;
    }
    mean /= norm;
    for (size_t i = 0; i < input_->size(); ++i) {
      for (int tr = 0; tr < tracks_.size(); ++tr) {
        v(tr) = (double) tracks_[tr]->get(i);

      }
      stddev += weights(k,i) * (v - mean) * (v - mean).transpose() ;
    }
    stddev /= norm;    
    emits_[k].set_mean(mean);
    emits_[k].set_var(stddev);
    //    std::cerr << "State " << k << std::endl;
    //    std::cerr << "Mean " << mean  << std::endl;
    //    std::cerr << "Std " << stddev << std::endl;
  }        
}

void 
MVGaussMultiTrackHMM::UpdateEmissionDistGibbs(Rng& r, const StateVectorType& states)
{
  ERRORLOG("NOT IMPLMEENTED YET");
}

void
MVGaussMultiTrackHMM::UpdateSoftEvidence(MatrixType& softev)
{  
  if (softev.rows() != num_states_ || softev.cols() != (int)input_->size())
    softev.resize(num_states_, input_->size());
  int nt = tracks_.size();
  Eigen::VectorXd v(nt);
  for (int i = 0; i < softev.cols(); ++i) {
    for (int j = 0; j < softev.rows(); ++j) {
      for (int k = 0; k < nt; ++k) {
        v(k) = (double)tracks_[k]->get(i);
      }
      softev(j, i) = emits_[j].Pdf(v) + 1.11e-16;
    }
  }
}



