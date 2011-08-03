#include "analysis/Dist.h"

double
GaussDist::Pdf(double val)
{
  return (gsl_ran_gaussian_pdf(val - m_, stddev_));
}

double
GaussDist::Sample(const gsl_rng* r)
{
  return (m_ + gsl_ran_gaussian(r, stddev_));
}

std::vector<double>
GaussDist::Sample(int n, const gsl_rng* r)
{
  std::vector<double> v(n);
  for (int i = 0; i < n; ++i) {
    v[i] = Sample(r);
  }
  return v;
}


//----------------------------------------------------------

double 
BernDist::Pdf(int k)
{
  return gsl_ran_bernoulli_pdf(k, p_);
}

int 
BernDist::Sample(const gsl_rng* r)
{
  return gsl_ran_bernoulli(r, p_);
}

//----------------------------------------------------------



double
MVGaussDist::Pdf(const Eigen::VectorXd& pt)
{
  double quadform = (m_ - pt).transpose() * varI_ * (m_ - pt);
  double logl = -0.5 * quadform - log(sqrt(vardet_)) - ndim_*log(2*M_PI)/2.0;
  return exp(logl);
}

void
MVGaussDist::set_var(const Eigen::MatrixXd& var) 
{
  assert(var.rows() == ndim_ && var.cols() == ndim_);
  // var = U'U, where U is upper triangular
  Eigen::MatrixXd cholfact_ = var.llt().matrixU(); 
  Eigen::MatrixXd ttI = cholfact_.inverse();
  varI_ = ttI * ttI.transpose();
  vardet_ = 1.0;
  // p 147 of MLaPP
  // |Sigma| = \prod_{i=1}^D R_{i,i}^2
  for (int i = 0; i < var.cols(); ++i) {
    vardet_ *= cholfact_(i,i)*cholfact_(i,i);
  }
}

Eigen::VectorXd
MVGaussDist::Sample(gsl_rng* r)
{
  GaussDist g(0.0, 1.0);
  Eigen::VectorXd temp(ndim_);
  for (int i = 0; i < ndim_; ++i) {
    temp(i) = g.Sample(r);
  }
  return (m_ + cholfact_.transpose() * temp);
}
       
//----------------------------------------------------------


//----------------------------------------------------------

//----------------------------------------------------------


MultiDist::MultiDist()
: n_(1)
, lookup(NULL)
{
  set_vals(Eigen::VectorXd::Constant(n_, 1.0/((double) n_)));
}

MultiDist::MultiDist(int n)
: n_(n)
, lookup(NULL)
{

  set_vals(Eigen::VectorXd::Constant(n_, 1.0/((double) n)));
}

void
MultiDist::set_vals(const Eigen::VectorXd& vals)
{
  if (vals.size() != n_) {
    n_ = vals.size();
  }
  if (lookup != NULL)
    gsl_ran_discrete_free(lookup);
  double* p = new double[n_];
  for (int i = 0; i < n_; ++i) {
    p[i] = vals(i);
  }
  lookup = gsl_ran_discrete_preproc(n_, p);
  delete[] p;
}

double
MultiDist::Pdf(int val)
{
  return gsl_ran_discrete_pdf((size_t)val, lookup);
}

int
MultiDist::Sample(const gsl_rng* r)
{
  return ((int) gsl_ran_discrete(r, lookup));
}

//----------------------------------------------------------

double
ChiSqDist::Sample(gsl_rng* r)
{
  return gsl_ran_chisq(r, nu_);
}

double
ChiSqDist::Pdf(double x)
{
  return gsl_ran_chisq_pdf(x, nu_);
}

//----------------------------------------------------------

double
StudentTDist::Sample(gsl_rng* r)
{
  // \mu + \sigma * z * \sqrt{\nu / x}
  double z = gauss_.Sample(r);
  double x = chi_.Sample(r);
  return (gauss_.mean() + sqrt(gauss_.stddev()) * z * sqrt(chi_.df() / x));
}


//----------------------------------------------------------

double
InvChiSqDist::Sample(gsl_rng* r)
{
  return (base_dist_.df() * s2_) / base_dist_.Sample(r);
}


//----------------------------------------------------------

DirichletDist::DirichletDist() : K_(1)
{
  alpha_ = new double[1];
}

DirichletDist::DirichletDist(int K) : K_(K)
{
  alpha_ = new double[K];
  for (int i = 0; i < K; ++i) {
    alpha_[i] = 1.0/((double) K);
  }
}

void
DirichletDist::set_alpha(const Eigen::VectorXd& alpha)
{
  if (K_ != alpha.size()) {
    K_ = alpha.size();
    delete[] alpha_;
    alpha_ = new double[K_];
  }
  for (int i = 0; i < K_; ++i) {
    alpha_[i] = alpha(i);
  }
}

Eigen::VectorXd
DirichletDist::alpha()
{
  Eigen::VectorXd v(K_);
  for (int i = 0; i < K_; ++i) {
    v(i) = alpha_[i];
  }
  return v;
}

Eigen::VectorXd 
DirichletDist::Sample(gsl_rng* r)
{
  double* temp = new double[K_];
  Eigen::VectorXd ret(K_);
  gsl_ran_dirichlet(r, (size_t)K_, alpha_, temp);
  for (int i = 0; i < K_; ++i) {
    ret(i) = temp[i];
  }
  delete[] temp;
  return ret;
}

double 
DirichletDist::Pdf(const Eigen::VectorXd& x)
{
  assert(x.size() == K_);
  const double* d = x.data();
  return gsl_ran_dirichlet_pdf(K_, alpha_, d);
}
