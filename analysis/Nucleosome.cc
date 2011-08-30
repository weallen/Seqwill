#include "analysis/Nucleosome.h"

void
NucPosFinder::InitializeEvents()
{
    
}


// each read maintains a distribution over 
// every possible location
void
NucPosFinder::FitEM() 
{
    // E STEP
    float ll = 0.0; 
    for (size_t i = 0; i < frags_.size(); ++i) {
        
    }
    
    float prior_ll = 0.0;
    for (size_t i = 0; i < events_.size(); ++i) {
        prior_ll += events_[i].prob;
    }
    prior_ll *= alpha_;
    ll -= prior_ll;
    
    // M STEP
    
}

void
NucPosFinder::UpdateEmpiricalDist()
{
    Eigen::VectorXd dist = Eigen::VectorXd::Zero(width_);
    
    for(std::vector<Event>::iterator it = events_.begin();
        it != events_.end(); ++it) {
        int start = (int)std::max(0, it->pos - width_);
        int stop = (int)std::min(pileup_->stop(), it->pos + width_);
        for (int i = start; i < stop; ++i) {
            int temp = pileup_->get(i);
            dist(abs(it->pos - i)) += temp;
        }
    }
    
    dist /= dist.sum();
    dist_.set_vals(dist);
}

//-------------------------------------------------

void
NucPileup::ComputeProcess()
{    
    output_->set_trackname(tname_);
    output_->set_subtrackname(stname_);
    output_->set_extends(start_, stop_);
    for (size_t i = 0; i < output_->size(); ++i) {
        output_->set(i, 0);
    }
    int nfrags = 0;
    float mean_fraglen = 0.0;
    float var_fraglen = 0.0;

    // compute mean and var of distribution over fragment lengths
    for (SingleReadFactory::iterator it = reads_->begin();
         it != reads_->end(); ++it) {
        if (it->second.partner_ref == it->second.ref_id) {
            nfrags++;
            mean_fraglen += (float)abs(it->second.pos - it->second.partner_pos);
            int pos = (int)floor((float)(it->second.partner_pos + it->second.pos)/2.0);
            int temp = output_->get(pos);
            output_->set(pos, temp+1);
        }        
    }
    mean_fraglen /= nfrags;

    for (SingleReadFactory::iterator it = reads_->begin();
         it != reads_->end(); ++it) {
        if (it->second.partner_ref == it->second.ref_id) {
            float len = (float)abs(it->second.pos - it->second.partner_pos);
            var_fraglen += pow(len - mean_fraglen,2);
        }                
    }
    var_fraglen /= nfrags;
    mean_ = mean_fraglen;
    var_ = var_fraglen;
}
 


//-------------------------------------------------

void
NucKDE::ComputeAnalysis()
{
    output_->set_extends(0, input_->size());
    output_->set_trackname(tname_);
    output_->set_subtrackname(stname_);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, input_->size()),
                      TBB_NucKernel(*input_, *output_, w_));
    
}

//-------------------------------------------------

void
TBB_NucKernel::operator()(const tbb::blocked_range<size_t>& r) const
{
    for (size_t pos = r.begin(); pos != r.end(); ++pos) {
        float val = D(pos);
        
        int start = std::max(0, (int)pos - 150);
        int stop = std::min(track_.size(), pos + 150);
        float denom = 0.0;
        for (int i = start; i < stop; ++i) {
            denom +=  D(i);
        }
        denom *= (1.09 / w_) * (stop - start);        
        output_.set(pos, val / denom); 
    }
}



float 
TBB_NucKernel::D(int i) const
{
    float val = 0.0;
    int start = std::max(0, (int)i - w_);
    int stop = std::min((int)i + w_, track_.stop());
    for (int j = start; j < stop; ++j) {
        float temp = (float) i - j;
        val += pow(1.0 - (temp / (float)w_)*(temp / (float)w_), 3);
        val *= (float)track_[j];
    }
    return val;
}