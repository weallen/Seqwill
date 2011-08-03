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
MVGaussDist::Pdf(const VectorType& pt)
{
  double denom = sqrt(pow(2 * M_PI, ndim_) * var_.determinant());
  double temp = (pt - m_).transpose() * var_.inverse() * (pt - m_);
  double numer = exp(-(0.5 * temp));
  return (numer / denom);
}

MVGaussDist::VectorType
MVGaussDist::Sample(gsl_rng* r)
{
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
DirichletDist::set_alpha(const VectorType& alpha)
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

DirichletDist::VectorType
DirichletDist::alpha()
{
  VectorType v(K_);
  for (int i = 0; i < K_; ++i) {
    v(i) = alpha_[i];
  }
  return v;
}

DirichletDist::VectorType 
DirichletDist::Sample(gsl_rng* r)
{
  double* temp = new double[K_];
  VectorType ret(K_);
  gsl_ran_dirichlet(r, (size_t)K_, alpha_, temp);
  for (int i = 0; i < K_; ++i) {
    ret(i) = temp[i];
  }
  delete[] temp;
  return ret;
}

double 
DirichletDist::Pdf(const VectorType& x)
{
  assert(x.size() == K_);
  const double* d = x.data();
  return gsl_ran_dirichlet_pdf(K_, alpha_, d);
}
