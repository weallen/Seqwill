#ifndef MEDIP_NORMALIZE_H_
#define MEDIP_NORMALIZE_H_
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <math.h>

#include <Eigen/Dense>

#include "base/MatrixUtil.h"

#include "analysis/AnalysisBase.h"
#include "common/Track.h"
#include "analysis/Dist.h"

/*
 * 1. set all cpgs to 1
 * repeat:
 *   2. assign fractions of fragments to CpGs based on scores
 *   3. sum fractions of reads for each CpG
 * for each fragment
 */


class MedipNormalize : public Analysis<float, float>
{
public:
  using Analysis<float, float>::analysis_name_;
  
  MedipNormalize()    
  : frag_len_(-1) 
  , alpha_(std::numeric_limits<double>::quiet_NaN())
  , beta_(std::numeric_limits<double>::quiet_NaN())
  { analysis_name_ = std::string("MedipNormalize"); }

  virtual ~MedipNormalize() {}

  void set_frag_len(int frag_len)
  { frag_len_ = frag_len; }
  
  void set_chr(Track<unsigned char>::Ptr chr)
  { chr_ = chr; }
  
  const Eigen::VectorXd& coupling() const
  { return coupling_; }
  
private:  
  void ComputeCoupling();
  virtual void ComputeAnalysis();
  void Sample();
  void FitBeta(const std::vector<double>& calls, BetaDist& b);
  void FitLinear();
  double SampleMethState(double methval, double coupling);
                         
  Track<unsigned char>::Ptr chr_;
  int frag_len_;
  // Phi(i) = alpha_ + beta_*X(i)
  double alpha_;
  double beta_;

  Eigen::VectorXd coupling_;
};



// divides genome into intervals and counts CpGs 
// (possibly and CpAs) in each interval
class CpGCounter : public Analysis<unsigned char, int>
{
public:
  CpGCounter() 
  : do_cpa_(false)
    , res_(1)
    , tname_("")
    , stname_("")
  { analysis_name_ = std::string("CpGCounter"); }
  
  virtual ~CpGCounter() {}
  
  void set_resolution(int res) 
  { assert(res > 0); res_ = res; }
  
  void set_do_cpa(bool t) 
  { do_cpa_ = t; }
 
  
private:
  virtual void ComputeAnalysis();
  bool do_cpa_;
  int res_;
  std::string tname_;
  std::string stname_;
};
#endif
