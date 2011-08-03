#include "analysis/Kmeans.h"

void
Kmeans::add_track(Track<float>::Ptr track)
{
  num_tracks_++;
  tracks_.push_back(track);
  Init();
}

void
Kmeans::Init() 
{
  for (int k = 0; k < K_; ++k) {    
    means_[k] = Eigen::VectorXd::Random(num_tracks_);
    vars_[k] = Eigen::MatrixXd::Identity(num_tracks_, num_tracks_);
  }
}

std::vector<Eigen::VectorXd>
Kmeans::Fit() 
{

  double delta = INFINITY;
  Eigen::VectorXd curr(num_tracks_);
  Eigen::VectorXi assignments(tracks_[0]->size());
  Eigen::VectorXd dists(K_);
  int curr_assign;
  Eigen::VectorXi curr_assign_count;
  Eigen::MatrixXd curr_assign_sum;

  int n = 0;
  while (n < 10) {
    n++;
    curr_assign_count = Eigen::VectorXi::Zero(K_);
    curr_assign_sum = Eigen::MatrixXd::Zero(num_tracks_, K_);
    DEBUGLOG("Kmeans step " + Stringify(n));
    for (int i = 0; i < (int)tracks_[0]->size(); ++i) {
      for (int t = 0; t < num_tracks_; ++t) {
	curr(t) = tracks_[t]->get(i);
      }
      for (int k = 0; k < K_; ++k) {
	dists(k) = (curr - means_[k]).norm();
      }
      dists.minCoeff(&curr_assign);
      assignments(i) = curr_assign;
    }

    for (int k = 0; k < K_; ++k) {
      for (int t = 0; t < (int)tracks_[0]->size(); ++t) {
	curr_assign_count(assignments(t)) += 1;
	for (int i = 0; i < num_tracks_; ++i) {
	  curr_assign_sum(i, k) += tracks_[i]->get(t);
	}
      }
      means_[k] =  curr_assign_sum.col(k) / curr_assign_count(k);
    }
  }

  return means_;
}
