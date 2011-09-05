#include "analysis/Kmeans.h"

void
Kmeans::set_track(Track<float>::Ptr track)
{
  track_ = track;
  Init();
}

void
Kmeans::Init() 
{
  Rng* r = new Rng;
  for (int k = 0; k < K_; ++k) {    
    means_[k] = gsl_rng_uniform(r->rng());
  }
  delete r;
}

void
Kmeans::Fit() 
{
  double* mean_sums = new double[K_];
  int* mean_counts = new int[K_];
  double* old_means = new double[K_];

    // initialize with random paritions
  for (size_t i = 0; i < track_->size(); ++i) {
    int idx = rand() % K_;
    mean_sums[idx] += track_->get(i);
    mean_counts[i] += 1;
  }

  for (int i = 0; i < K_; ++i) {
      means_[i] = mean_sums[i] / (double)mean_counts[i];
      old_means[i] = 0.0;
  }

  // fit model
  double delta = INFINITY;
  int n = 0;
  while (delta > 1e-4) {
    n++;
    for (int i = 0; i < K_; ++i) {
      mean_sums[i] = 0.0;
      mean_counts[i] = 1;
    }
    for (size_t i = 0; i < track_->size(); ++i) {
      double min_dist = 1e50;
      int min_idx = 0;
      double temp_dist;
        double curr = (double)track_->get(i);
      for (int k = 0; k < K_; ++k) {
        temp_dist = pow(curr - means_[k],2);
        if (temp_dist < min_dist) {
          min_dist = temp_dist;
          min_idx = k;
        }
      }
      mean_sums[min_idx] += curr;
      mean_counts[min_idx] += 1;
    }

    delta = 0.0;
    for (int i = 0; i < K_; ++i) {
      means_[i] = mean_sums[i] / (double)mean_counts[i];
      delta += abs(means_[i] - old_means[i]);
        old_means[i] = means_[i];
    }
  }
  delete old_means;
  delete mean_sums;
  delete mean_counts;
}
