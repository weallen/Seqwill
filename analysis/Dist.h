#ifndef HMM_EMIT_H_
#define HMM_EMIT_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <gsl/gsl_vector.h>

#include <math.h>
#include <cassert>
#include "analysis/gsl_addon.h"

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

class MVGaussDist
{
public:  
  typedef Eigen::VectorXf VectorType;
  typedef Eigen::MatrixXf MatrixType;

  MVGaussDist()
  : ndim_(2)
  , m_(2)
  , var_(2,2)
  {  }
  
  MVGaussDist(const VectorType& m, const MatrixType& c)
  : ndim_(m.size())
  {
    set_mean(m);
    set_var(c);
  }
  virtual ~MVGaussDist(){}
  
  void set_ndim(int n)
  { ndim_ = n; m_.resize(n); var_.resize(n,n); }
  
  const int ndim() const
  { return ndim_; }
  
  void set_mean(const VectorType& m)
  { m_ = m; }
  
  VectorType mean()
  { return m_; }
  
  void set_var(const MatrixType& c)
  { var_ = c; }
  
  MatrixType var()
  { return var_; }
  
  double Pdf(const VectorType& pt);
  VectorType Sample(gsl_rng* r); 

private:
  int ndim_;
  VectorType m_;
  MatrixType var_;
};

class MVStudentTDist
{
public:
  typedef Eigen::VectorXf VectorType;
  typedef Eigen::MatrixXf MatrixType;

  MVStudentTDist();
  virtual ~MVStudentTDist();
  
  void set_loc(const VectorType& loc);
  VectorType loc() const;
  void set_scale(const MatrixType& scale);
  MatrixType scale() const;
  void set_dof(int dof);
  const int dof() const;
  
  double Pdf(VectorType pt);
  VectorType Sample(gsl_rng* r);

private:
  gsl_vector* location_;
  gsl_matrix* scale_;
  int dof_;
};

// Just for sampling for now
class WishartDist
{
public:
  typedef Eigen::VectorXf VectorType;
  typedef Eigen::MatrixXf MatrixType;

  WishartDist();
  virtual ~WishartDist();
  
  double Pdf(const VectorType& pt);
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
  
  void set_scale(double s2) { s2_ = s2; }
  void set_df(int n) { base_dist_.set_df(n); }
  
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

