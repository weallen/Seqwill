#include "analysis/HMM.h"

float 
HMM::LogProb()
{
  float logprob = 0.0;
  for (int i = 0; i < emit_probs.size(); ++i) {
    logprob += log(soft_evidence_.col(i).maxCoeff());
  }
  return logprob;
}

void 
HMM::FitEM() 
{
}

void 
HMM::FwdBack()
{
  
}

void
HMM::ViterbiDecode(StateVectorType& path)
{ 
  int T = ((int) input_->size()) - 1;
  
  // Delta holds the normalized values
  MatrixType delta = MatrixType::Zero(num_states_, T);
  
  // psi stores the indices of the maximum value for each
  StateMatrixType psi = MatrixType::Zero(num_states_, T);
  
  // Initialize last element as soft evidence weighted
  // By prior probabilities
  delta.col(T) = init_ * soft_evidence_.col(T);  
  VectorType v;
  MatrixType::Index idx;
  for (int t=1; t <= T; ++t) {
    for (int j=0; j < num_states_; ++j) {
      // find the unnormalized values
      v = delta.col(t-1) * trans_.col(j);
      delta(j,t) = v.maxCoeff(&idx) * soft_evidence_(j,t);
      psi(j,t) = (int) idx;
    }
    // Normalize 
    delta.col(t) /= delta.col(t).sum();
  }
  
  // Traceback
  // First find the 

  float maxcoeff = delta.col(T).maxCoeff(&idx);
  path(t) = idx;  
  for (t = T - 1; t >= 0; --t) {
    path(t) = psi(path(t+1), t+1);
  }
}


HMM::MatrixType 
HMM::CountTransitions(const StateVectorType& curr_states )
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

//----------------------------------------------------------------