#ifndef KMEANS_H_
#define KMEANS_H_
#include <iostream>
#include <Eigen/Dense>
#include "common/Track.h"
#include "base/StringUtil.h"
#include "analysis/AnalysisBase.h"

class Kmeans
{
public:
  Kmeans() : K_(1), num_tracks_(0), means_(K_), vars_(K_) {  }
  Kmeans(int k) : K_(k), num_tracks_(0), means_(K_), vars_(K_) { }
  virtual ~Kmeans() {}
  
  void add_track(Track<float>::Ptr track);

  std::vector<Eigen::VectorXd> Fit();

private:
  void Init(); 

  int K_;
  int num_tracks_;
  std::vector<Eigen::VectorXd> means_;
  std::vector<Eigen::MatrixXd> vars_; 
  std::vector<Track<float>::Ptr> tracks_;
};

#endif
