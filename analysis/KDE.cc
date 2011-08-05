#include "analysis/KDE.h"

float
KDE::Eval(float x)
{
  float output = 0.0;
  for (int i = 0; i < data_.size(); ++i) {
    output += kernfn_((x - data_(i)), bandwidth_);
  }
  output /= (data_.size() * bandwidth_);
  return output;
}

float
GaussianKernel(float x, float b)
{
  float output = exp(-(x*x) * 0.5) / sqrt(2*M_PI);
  return (output / b);
}

