#ifndef KMEANS_H_
#define KMEANS_H_
#include <iostream>
#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <math.h>

#include "common/Track.h"
#include "base/StringUtil.h"
#include "analysis/AnalysisBase.h"
#include "analysis/Random.h"
class Kmeans
{
public:
  Kmeans() : K_(1), means_(K_) {}
  Kmeans(int k) : K_(k), means_(K_) {}
  virtual ~Kmeans() { }
  
  void set_track(Track<float>::Ptr track);

  void Fit();

  std::vector<double> means()
  { return means_; }

private:
  void Init(); 

  int K_;
  std::vector<double> means_;
  Track<float>::Ptr track_;
};

#endif
