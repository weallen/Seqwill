#ifndef RANDOM_H_
#define RANDOM_H_
#include <gsl/gsl_rng.h>

class Rng 
{
public:
  Rng() { Init(); }
  virtual ~Rng() { gsl_rng_free(rng_); }
  
  void Init();
  
  gsl_rng* rng() 
  { return rng_; }
    
private:
  gsl_rng* rng_;
};

#endif
