#ifndef HMM_EMIT_H_
#define HMM_EMIT_H_

#include <vector>
#include <Eigen/Dense>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <cassert>

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

  double mean() const { return m_; }

  void set_stddev(float s) 
  { stddev_ = s; }

  float stddev() const
  { return stddev_; }

  double Pdf(double val);  
  double Sample(const gsl_rng* r);
  std::vector<double> Sample(int n, const gsl_rng* r);

private:
  double m_;
  double stddev_;
};

class BernDist
{
public:
  BernDist() : p_(0.0) {}
  BernDist(double p) : p_(p) {}

  void set_prob(double p)
  { p_ = p; }
  
  const double prob() const
  { return p_; }
  
  double Pdf(int k);
  int Sample(const gsl_rng* r);
  
private:
  double p_;
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

  double pdf(VectorType pt) { return 0.0; }

  VectorType Sample() 
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
  MultiDist();
  MultiDist(int n);
  virtual ~MultiDist()
  { gsl_ran_discrete_free(lookup); }
 
  int num() const 
  { return n_; }

  void set_vals(const Eigen::VectorXd& vals);
  
  double Pdf(int val);
  
  int Sample(const gsl_rng* r);

private:
  int n_;
  gsl_ran_discrete_t* lookup;
};

// Used for estimating mean of gaussian
// for n independent Gaussian rvs with unit variance Y_i
// X_i = \sum_i Y_i^2 has X^2 distribution with n d.o.f.
class ChiSqDist
{
public:
  ChiSqDist() : nu_(1.0) {}
  ChiSqDist(double nu) : nu_(nu) {}
  virtual ~ChiSqDist() {}
  void set_df(double n) { nu_ = n; }
  const double df() const { return nu_; }
  double Sample(gsl_rng* r);
  double Pdf(double x);
  
private:
  double nu_;
};

// Just a sampling distribution for now...
class InvChiSqDist
{
public:
  InvChiSqDist() : base_dist_(1.0), s2_(1.0) {}
  InvChiSqDist(double nu, double scale) 
  : base_dist_(nu), s2_(scale)
  {}
  
  double scale() { return s2_; }
  double df() { return base_dist_.df(); }

  double Sample(gsl_rng* r);
  
private:
  ChiSqDist base_dist_;
  double s2_;
};

// Just a sampling distribution for now...
class StudentTDist
{
public:
  StudentTDist() {}
  StudentTDist(double nu, double mu, double sigma2) 
  : chi_(nu), gauss_(mu, sigma2) 
  {}
  
  virtual ~StudentTDist() {}
  
  void set_df(double n) { chi_.set_df(n); }
  void set_mean(double mean) { gauss_.set_mean(mean); }
  void set_stddev(double s) { gauss_.set_stddev(s); }
  const double df() const { return chi_.df(); }
  const double mean() const { return gauss_.mean(); }
  const double stddev() const { return gauss_.stddev(); }
  
  double Sample(gsl_rng* r);
  
private:
  ChiSqDist chi_;
  GaussDist gauss_;
};

class DirichletDist
{
public:
  typedef Eigen::VectorXd VectorType;
  
  DirichletDist();
  DirichletDist(int K);
  virtual ~DirichletDist() 
  { delete[] alpha_; }
  
  void set_alpha(const VectorType& alpha);
  VectorType alpha();
  
  VectorType Sample(gsl_rng* r);
  double Pdf(const VectorType& x);
  
private:
  double* alpha_;
  int K_;
  
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
 
private:
  int nmix_; // num of mixture components
  std::vector<VectorType> means_;
  std::vector<MatrixType> covs_;
};



#endif

