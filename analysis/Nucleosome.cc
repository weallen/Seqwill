#include "analysis/Nucleosome.h"



void
NucMapper::Compute()
{
    output_->set_extends(0, input_->size());
    output_->set_trackname(tname_);
    output_->set_subtrackname(stname_);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, input_->size()),
                      TBB_NucKernel(*input_, *output_));
    
}

//-------------------------------------------------

void
TBB_NucKernel::operator()(const tbb::blocked_range<size_t>& r) const
{
    for (size_t pos = r.begin(); pos != r.end(); ++pos) {
        float val = 0.0;
        for (size_t i = 0; i < track_.size(); ++i) {
                
        }
        output_.set(pos, val); 
    }
}

