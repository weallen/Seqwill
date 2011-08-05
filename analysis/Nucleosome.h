#ifndef NUCLEOSOME_H_
#define NUCLEOSOME_H_

#include <Eigen/Dense>
#include "io/BamIO.h"
#include "io/TrackIO.h"
#include "common/Track.h"

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
  
  const void num_bins() const
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

// stores one track for plus strand
// and one strand for minus strand for each chr
class NucData
{
public:
  
private:
};

class NucMapper
{
public:
  NucMapper() {}
  virtual ~NucMapper() {}
  
  
  void set_plus_track(Track<int>::Ptr plus_strand)
  { plus_ = plus_strand; }

  void set_minus_track(Track<int>::Ptr minus_strand)
  { minus_ = minus_strand; }

private:
  Track<int>::Ptr plus_;
  Track<int>::Ptr minus_;
};
#endif
