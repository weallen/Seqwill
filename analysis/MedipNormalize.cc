#include "analysis/MedipNormalize.h"


// 1. find maxima of CpG amounts
// Do linear regression only up to maxima
// Use linear regression values to sample true values
// for all

void
MedipNormalize::Calibrate()
{
  assert(frag_len_ != -1);
  
  // compute coupling factors
  coupling_ = Eigen::VectorXf::Zero(input_->size());
  int chrlen = chr_->stop();
  int res = input_->resolution();
  for (size_t i = 0; i < input_->size(); ++i) {
    int absmax = (i+1) * res;
    int absmin = i * res;
    int maxpos = std::min(absmax + frag_len_ - 1, chrlen-1);
    int minpos = std::max(absmin - frag_len_ + 1, 0);
    int total = 0;
    for (int j = minpos; j <= maxpos; ++j) {
      if (chr_->get(j) == 'C' && chr_->get(j+1) == 'G') {
	total++;
      }
    }
    coupling_(i) = (float)total;
  }
}

void
MedipNormalize::ComputeAnalysis()
{
  assert(!isnan(beta_) && !isnan(alpha_) && frag_len_ != -1);
  DEBUGLOG("Calibrating MedipNormalize...");
  Calibrate();
  DEBUGLOG("Fitting linear regression...");
  FitLinear();  
}
// A_p = number of reads overlapping a particular bin
// A_base = A_p-intercept of linreg
// r is slope of model
// sample from f(A|m) = \Pi_p G(A_p|A_base + r*C_cp, v^-1)
void
MedipNormalize::Sample()
{
}

void
MedipNormalize::FitLinear(float& alpha, float& beta)
{
  // find local maxima of 
  Eigen::VectorXf predict;
  Eigen::VectorXf response;

  // count number with coupling factors greater
  // than 10
  int count = 0;
  
  // Do LinReg
  SimpleLinReg(predict, response, beta);
  alpha = response.sum()/((float)response.size());
}

//------------------------------------------------------------------------------

void
FindCouplingFactors::ComputeAnalysis()
{
}

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
  unsigned char curr_char = '';
  unsigned char prev_char = '';
  for (size_t i = 0; i < output_->size(); ++i) {
    output_->set(i, 0);
  }
  
  int curr_count;
  for (size_t i = 0; i < input_->size(); i += res_) {
    curr_count = 0;
    prev_char = input_->get(i);
    for (size_t j = 1; j < res_; ++j) {
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
