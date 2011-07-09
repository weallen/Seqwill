#ifndef HMM_EMIT_H_
#define HMM_EMIT_H_

#include <Eigen/Dense>

template <int Ndim>
class MixGaussDist
{
public:
  typedef Eigen::Matrix<float, Ndim, 1> VectorType;
  typedef Eigen::Matrix<float, Ndim, Ndim> MatrixType;

  MixGaussDist() {}
  virtual ~MixGaussDist() {
  }

  void Init()
  {
    means_.resize(nmix_);
    covs_.resize(nmix_);
  }
 
  float operator()(float val)
  {}

private:
  int nmix_; // num of mixture components
  std::vector<VectorType> means_;
  std::vector<MatrixType> covs_;
};

template <int Ndim>
class GaussDist
{
public:
  typedef Eigen::Matrix<float, Ndim, 1> VectorType;

  GaussDist() {}
  virtual ~GaussDist() {}

  void Init();

  float operator()(float val)
  {}
  

};

template <int Ndim>
class MultiDist
{
public:
  MultiEmit() {}
  virtual ~MultiEmit() {}

  void Init(); 

  float operator()(int val)
  {}
};


#endif

