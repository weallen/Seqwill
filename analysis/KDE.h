#ifndef KDE_H_
#define KDE_H_

#include <Eigen/Dense>
#include <boost/function.hpp>
#include <math.h>

class KDE 
{
public:
  typedef boost::function2<float,float,float> KernFnType;

  KDE()
    : bandwidth_(1.0)
    , data_(1)
    {}

  KDE(const Eigen::VectorXf& v, 
      const KernFnType& f,
      float bandwidth) 
    : data_(v) 
    , kernfn_(f)
    , bandwidth_(bandwidth)
  { }

  virtual ~KDE() {}

  float Eval(float x);

  void set_data(const Eigen::VectorXf& v) 
  { data_ = v; }

  const Eigen::VectorXf& data() const
  { return data_; }
 
  void set_kern(const KernFnType k)
  { kernfn_ = k; }

  const KernFnType& kern() const 
  { return kernfn_; }
 
  void set_bandwidth(float b) 
  { bandwidth_ = b; }
  
  float bandwidth() 
  { return bandwidth_; }

private:
  Eigen::VectorXf data_;
  KernFnType kernfn_;
  float bandwidth_;
};

boost::function2<float,float, float> GaussKernFn;

float GaussianKernel(float x, float b);



#endif
