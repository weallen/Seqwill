#include "analysis/HSMM.h"




//-------------------------------------------------------------------

HSMM::~HSMM()
{
  delete duration_;
  delete obs_;
  delete trans_;
  delete init_;
  delete states_;
}



void
HSMM::FitEM()
{
    
}

void 
HSMM::FitBlockedGibbs()
{
    
}

void
HSMM::FilterFwd()
{
    
}