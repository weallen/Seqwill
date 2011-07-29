#ifndef HMM_EMIT_H_
#define HMM_EMIT_H_

#include <vector>
#include <Eigen/Dense>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

typedef Eigen::VectorXf VectorType;
typedef Eigen::MatrixXf MatrixType;

float SampleGauss(float m, float s);
int SampleMulti(Eigen::VectorXf vals);
VectorType SampleMultiVarGauss(VectorType m, MatrixType cov);


class GaussDist 
{
public:

  GaussDist() : m_(0.0), stddev_(1.0) {}
  GaussDist(float m, float s) : m_(m), stddev_(s) {}
  
  virtual ~GaussDist() {}

  void Init();

  void set_mean(float m)
  { m_ = m; }

  float mean() const { return m_; }

  void set_stddev(float s) 
  { stddev_ = s; }

  float stddev() const
  { return stddev_; }

  float pdf(float val);  
  float sample(const gsl_rng* r);
  std::vector<float> sample(int n, const gsl_rng* r);

private:
  float m_;
  float stddev_;
};

class MultiVarGaussDist
{
public:
  MultiVarGaussDist()
  : ndim_(0)
  {}
  virtual ~MultiVarGaussDist() {}
  
  void set_ndim(int n) 
  { ndim_ = n; }

  int ndim() const 
  { return ndim_; }

  void set_mean(const VectorType& m) 
  { m_ = m; }
  
  const VectorType& mean() const 
  { return m_; }

  void set_cov(const MatrixType& c) 
  { cov_ = c; }

  const MatrixType& cov() const
  { return cov_; }

  float pdf(VectorType pt) { return 0.0; }

  VectorType sample() 
  { VectorType v(10);  return v; }
  
  std::vector<VectorType> sample(int n);

private:
  int ndim_;
  VectorType m_;
  MatrixType cov_;
};

class MultiDist
{
public:
  typedef Eigen::VectorXf VectorType;

  MultiDist() : n_(1) {}
  virtual ~MultiDist() {}

  void Init(); 
 
  void set_num(int n) 
  { n_ = n; }

  int num() const 
  { return n_; }


  float pdf(int val)
  { return 0.0; }

  int sample() { return 0; }

private:
  int n_;
  VectorType vals_;
};

class MixGaussDist 
{
public:
  typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorType;
  typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> MatrixType;

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


#endif

