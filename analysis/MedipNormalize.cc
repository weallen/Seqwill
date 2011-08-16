#include "analysis/MedipNormalize.h"


// 1. find maxima of CpG amounts
// Do linear regression only up to maxima
// Use linear regression values to sample true values
// for all

void
MedipNormalize::ComputeCoupling()
{
  assert(frag_len_ != -1);
  
  // compute coupling factors
  coupling_ = Eigen::VectorXd::Zero(input_->size());
  int chrlen = chr_->stop();
  int res = input_->resolution();
  for (size_t i = 0; i < input_->size(); ++i) {
    int absmax = (i+1) * res;
    int absmin = i * res;
    int maxpos = std::min(absmax + frag_len_ - 1, chrlen-1);
    int minpos = std::max(absmin - frag_len_ + 1, 0);
    double total = 0.0;
    for (int j = minpos; j <= maxpos; ++j) {
      unsigned char curr = chr_->get(j);
      unsigned char next = chr_->get(j+1);
      if ((curr == 'C' && next == 'G')
          || (curr == 'c' && next == 'g')) {
        total++;
      }
    }
    coupling_(i) = total;
  }
}

void
MedipNormalize::ComputeAnalysis()
{
//  assert(!isnan(beta_) && !isnan(alpha_) && frag_len_ != -1);
  DEBUGLOG("Calibrating MedipNormalize...");
  ComputeCoupling();
  DEBUGLOG("Fitting linear regression...");
//  FitLinear();  
}

void
MedipNormalize::Sample()
{
  int N = (int)input_->size();
  double* array = new double[N];
  #pragma omp parallel for shared(array) private(i,j,b,calls)
  for (int i = 0; i < N; ++i) {
    std::vector<double> calls(100);
    for (int j = 0; j < 100; ++j) {
      calls[j] = SampleMethState(input_->get(i), coupling_(i));
    }
    BetaDist b;
    FitBeta(calls, b);
    array[i] = b.Mode();
  }  
}

// A_p = number of reads overlapping a particular bin
// A_base = A_p-intercept of linreg
// r is slope of model
// sample from f(A|m) = \Pi_p G(A_p|A_base + r*C_cp, v^-1)
//
// NESTED SAMPLING
// rewrite Z = \int p(D|\theta)p(\theta) d\theta
// as \int_0^1 L(\theta(x)) dx
// then use samples from prior with nested
// set of constraints of likelihood to evaluate that integral
// 1. sample elements from parameters space (sample from prior)
// 2. sort these elements by likelihood
// 3. 
double
MedipNormalize::SampleMethState(double methval, double coupling)
{
  
}

void 
MedipNormalize::FitBeta(const std::vector<double>& calls, BetaDist& b)
{
  
}

void
MedipNormalize::FitLinear()
{
  // find local maxima of 
  Eigen::VectorXd predict;
  Eigen::VectorXd response;

  // count number with coupling factors greater
  // than 10
  int count = 0;
  
  // Do LinReg
  SimpleLinReg(predict, response, beta_);
  alpha_ = response.sum()/((float)response.size());
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void 
CpGCounter::ComputeAnalysis() 
{
  assert(tname_ != "" && stname_ != "");
  output_ = TrackOutPtr(new Track<int>);
  output_->set_resolution(res_);
  output_->set_extends(floor(((float)input_->start()) / res_), floor(((float)input_->stop())/res_));
  output_->set_trackname(tname_);
  output_->set_subtrackname(stname_);
  unsigned char curr_char; 
  unsigned char prev_char;
  for (size_t i = 0; i < output_->size(); ++i) {
    output_->set(i, 0);
  }
  
  int curr_count;
  for (size_t i = 0; i < input_->size(); i += res_) {
    curr_count = 0;
    prev_char = input_->get(i);
    for (int j = 1; j < res_; ++j) {
      curr_char = input_->get(i+j);
      if (curr_char == 'G' && prev_char == 'C') 
	curr_count++;
      if (do_cpa_) {
	if (curr_char == 'A' && prev_char == 'C') {
	  curr_count++;
	}
      }
      prev_char = curr_char;
    }
    output_->set(i, curr_count);
  }
}
