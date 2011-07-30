#include "analysis/Dist.h"

double
GaussDist::pdf(double val)
{
  return (gsl_ran_gaussian_pdf(val - m_, stddev_));
}

double
GaussDist::sample(const gsl_rng* r)
{
  return (m_ + gsl_ran_gaussian(r, stddev_));
}

std::vector<double>
GaussDist::sample(int n, const gsl_rng* r)
{
  std::vector<double> v(n);
  for (int i = 0; i < n; ++i) {
    v[i] = sample(r);
  }
  return v;
}


//----------------------------------------------------------

double 
BernDist::pdf(int k)
{
  return gsl_ran_bernoulli_pdf(k, p_);
}

int 
BernDist::sample(const gsl_rng* r)
{
  return gsl_ran_bernoulli(r, p_);
}
