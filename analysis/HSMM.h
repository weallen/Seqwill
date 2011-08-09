#ifndef HSMM_H_
#define HSMM_H_

#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>

class HSMM
{
public:
  HSMM() {}
  virtual ~HSMM() {}
  void Resample();
  void Generate();
  
  void set_duration(Duration::Ptr d)
  { duration_ = d; }

  void set_trans(Transition::Ptr t)
  { }


  void Generate(Eigen::ArrayXi& stateseq);

private:
  void Resample(const Eigen::ArrayXXf obs);

  Duration::Ptr duration_;
  Obs::Ptr obs_;
  Transition::Ptr trans_;
  InitStates::Ptr init_;
  States::Ptr states_;
};

class InitStates
{
public:
  InitStates() {}
  virtual ~InitStates() {}
  
  void Resample(const Eigen::VectorXf& states);

private:
  int nstates_;  
  Eigen::VectorXf pi_;
};

class States
{
public:
  typedef boost::shared_ptr<States> Ptr;

  States() {}
  virtual ~States() {}
  
  void Resample();
  void Generate();
  void GenerateObs();
  void GenerateStates();

};

class Transition
{
public:
  typedef boost::shared_ptr<Transition> Ptr;

  Transition() {}
  virtual ~Transition() {}

private:
};

class Obs
{
public:
  typedef boost::shared_ptr<Obs> Ptr;

  Obs() {}
  virtual ~Obs() {}
  virtual void Resample(const Eigen::VectorXi& stateseq);
  
};

class GaussObs : public Obs
{
public:
  typedef boost::shared_ptr<GaussObs> Ptr;

  GaussObs() {}
  virtual ~GaussObs() {}
  virtual void Resample(const Eigen::VectorXi& stateseq);
};

class MVGaussObs : public Obs
{
public:
  typedef boost::shared_ptr<MVGaussObs> Ptr;

  MVGaussObs() {}
  virtual ~MVGaussObs() {}
  virtual void Resample(const Eigen::VectorXi& stateseq)
};

class Duration
{
public:
  typedef boost::shared_ptr<Duration> Ptr;

  Duration() {}
  virtual ~Duration() {}

};

class PoissonDuration : public Duration
{
public:
  typedef boost::shared_ptr<PoissonDuration> Ptr;

  PoissonDuration() {}
  virtual ~PoissonDuration() {}

private:
};


#endif
