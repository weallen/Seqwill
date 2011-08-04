#ifndef NUCLEOSOME_H_
#define NUCLEOSOME_H_

#include "common/Track.h"

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
