#ifndef NUCLEOSOME_H_
#define NUCLEOSOME_H_

#include <tbb/parallel_for.h>

#include <Eigen/Dense>
#include "io/BamIO.h"
#include "io/TrackIO.h"
#include "common/Track.h"
#include "base/Types.h"
#include "analysis/AnalysisBase.h"

class Histogram
{
public:
    // defaults to 10 bins
    Histogram() 
    : n_(10) 
    , counts_(10) {}
    
    Histogram(int n) 
    : n_(n)
    , counts_(n) {}
    
    virtual ~Histogram() {}
    
    void set_num_bins(int n)
    { n_ = n; }
    
    const int num_bins() const
    { return n_; }
    
    void add_to_bin(int i, int val) 
    { counts_[i] += val; }
    
    void set_bin(int i, int val) 
    { counts_[i] = val; }
    
    const int get_bin(int i) const
    { return counts_[i]; }
    
private:
    int n_; // num bins
    Eigen::ArrayXi counts_;
};


class NucMapper : public Analysis<PlusMinusDataInt, float>
{
public:
    NucMapper() {}
    virtual ~NucMapper() {}
    
    virtual void Compute();

private:
    
};

class TBB_NucKernel
{
public: 
    TBB_NucKernel(Track<PlusMinusDataInt>& track,
                  Track<float>& output)
    : track_(track)
    , output_(output)
    { }
    
    void operator()(const tbb::blocked_range<size_t>& r) const; 

private:
    Track<PlusMinusDataInt>& track_;
    Track<float>& output_;
};

#endif
