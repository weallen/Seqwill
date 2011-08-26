#ifndef HSMM_H_
#define HSMM_H_

#include <vector>

#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>

class Transition
{
public:
    Transition() {}
    virtual ~Transition() {}
    
    void set_trans(const Eigen::MatrixXd& trans)
    { trans_ = trans; }
    
    const Eigen::MatrixXd trans() const
    { return trans_; }
    
private:
    Eigen::MatrixXd trans_;
};

class Obs
{
public:
    Obs() {}
    
    virtual ~Obs() {}
    
    // XXX Not implemented yet
    virtual void Resample(const Eigen::VectorXi& stateseq) {}
    
};

class GaussObs : public Obs
{
public:
    GaussObs() {}
    virtual ~GaussObs() {}
    
    // XXX Not implemented yet
    virtual void Resample(const Eigen::VectorXi& stateseq) {}
};

class MVGaussObs : public Obs
{
public:
    MVGaussObs() {}
    virtual ~MVGaussObs() {}
    
    // XXX Not implemented yet
    virtual void Resample(const Eigen::VectorXi& stateseq) {}
};

class Duration
{
public:
    Duration() {}
    virtual ~Duration() {}
    
};

class GeometricDuration : public Duration
{
public:
    GeometricDuration() {}
    virtual ~GeometricDuration() {}
};

class PoissonDuration : public Duration
{
public:
    
    PoissonDuration() {}
    virtual ~PoissonDuration() {}
    
private:
    
};

class InitStates
{
public:
    InitStates() {}
    virtual ~InitStates() {}
    
    // XXX Not implemented yet
    void Resample(const Eigen::VectorXd& states) {}
    
private:
    int nstates_;  
    Eigen::VectorXd pi_;
};

class States
{
public:
    States() {}
    virtual ~States() {}
    

    /*
    void Resample();
    void Generate();
    void GenerateObs();
    void GenerateStates();
     */
private:    
    std::vector<double> durations_;
    std::vector<int> state_seq_;
};


class HSMM
{
public:
    HSMM() {}
    virtual ~HSMM();
    
    void Resample();
    void Generate();
    
    
    void FitEM();
    void FitBlockedGibbs();
    
    void FilterFwd(); 
    void SmoothBack();
    
    void set_duration(const Duration& d)
    { *duration_ = d; }
    
    void set_trans(const Transition& t)
    { trans_->set_trans(t.trans()); }

private:
    
    virtual void Generate(Eigen::ArrayXi& stateseq);
    virtual void Resample(const Eigen::ArrayXXd& obs);
    
    void Init();
    
    Duration* duration_;
    Obs* obs_;
    Transition* trans_;
    InitStates* init_;
    States* states_;
};


#endif
