#include "analysis/Dist.h"

float
GaussDist::pdf(float val)
{
  return ((float) gsl_ran_gaussian_pdf(val - m_, stddev_));
}

float
GaussDist::sample(const gsl_rng* r)
{
  return (m_ + (float) gsl_ran_gaussian(r, stddev_));
}

std::vector<float>
GaussDist::sample(int n, const gsl_rng* r)
{
  std::vector<float> v(n);
  for (int i = 0; i < n; ++i) {
    v[i] = sample(r);
  }
  return v;
}