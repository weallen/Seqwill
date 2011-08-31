#include "analysis/Nucleosome.h"

void
NucPosFinder::InitializeEvents()
{
    float num_events = 0.0;
    for (size_t i = 0; i < pileup_->size(); ++i) {
        if (pileup_->get(i) > 1) {
            Event* e = new Event;
            e->pos = (int)i;
            e->prob = 1.0;
            events_.push_back(e);
            num_events++;
        }
    }
    for (std::list<Event*>::iterator it = events_.begin();
         it != events_.end(); ++it) {
        (*it)->prob = 1.0/num_events;
    }
}

void 
NucPosFinder::AssignEventsToFrags()
{
    int num_bins = ceil((float) pileup_->size() / bin_size_);
    std::vector<std::vector<Event*> > idx(num_bins);
    for (std::list<Event*>::iterator it = events_.begin();
         it != events_.end(); ++it) {
        Event* event = *it;
        idx[PosToBinIdx(event->pos)].push_back(event);
    }
    for (std::vector<NucFrag>::iterator frag = frags_.begin();
         frag != frags_.end(); ++frag) {
        int start = std::max(PosToBinIdx(frag->start)-1,0);
        int stop = std::min(PosToBinIdx(frag->stop)+1, num_bins);
        for (int i = start; i <= stop; ++i) {
            if (idx[i].size() > 0) {
                std::vector<Event*> e = idx[i];
                for (std::vector<Event*>::iterator it = e.begin();
                     it != e.end(); ++it) {
                    Event* event = *it;
                    if (event->pos >= frag->start && event->pos <= frag->stop)
                        frag->nucs.push_back(event);
                }
            }
        }
        float num_nucs = (float) frag->nucs.size();
        for (std::vector<Event*>::iterator event = frag->nucs.begin();
             event != frag->nucs.end(); ++event) {
            frag->nuc_probs.push_back(1.0 / num_nucs);
        }
    }
}

void
NucPosFinder::ReadsToFrags()
{
    assert(reads_ != NULL);
    
    int nfrags = 0;
    // index frags by position of first mate
    for (SingleReadFactory::iterator it = reads_->begin();
         it != reads_->end(); ++it) {
        NucFrag f;
        //XXX may need to fix this
        if (it->second.partner_ref == it->second.ref_id) {
            if (it->second.pos > it->second.partner_pos) {
                f.stop = it->second.pos;
                f.start = it->second.partner_pos;
            } else {
                f.stop = it->second.partner_pos;
                f.start = it->second.pos;
            }          
        }
        frags_.push_back(f);
        nfrags++;
    }
    num_frags_ = nfrags;
}

// XXX NOT FULLY IMPLEMENTED YET
void
NucPosFinder::FitEM() 
{
    float delta_loglik = INFINITY;
    float old_loglik = 0.0;
    int n = 0;
    
    while (delta_loglik > 0.01) {
        DEBUGLOG("EM STEP " + Stringify(n));
        
        n++;
        
        // Init: Zero all the weights
        for (std::list<Event*>::iterator it = events_.begin(); 
             it != events_.end(); ++it) {
            (*it)->weight = 0.0;
        }
        
        // E STEP
        // update nucleosome positions
        
        
        // Estimate weights of events
        for (std::vector<NucFrag>::iterator frag = frags_.begin();
             frag != frags_.end(); ++frag) {
            for (size_t i = 0; i < frag->nucs.size(); ++i) {
                Event* event = frag->nucs[i];
                event->weight += frag->nuc_probs[i];
            }
        }
        
        // M STEP
        float ll = 0.0; 
        for (size_t i = 0; i < frags_.size(); ++i) {
            
        }
        
        // Update event probabilities
        float event_sum = 0.0;
        for (std::list<Event*>::iterator it = events_.begin();
             it != events_.end(); ++it) {
            Event* event = *it;
            event->prob = std::max(0.0, event->weight - alpha_);
            if (event->prob == 0.0) {
                events_.erase(it);
                delete event;
            } else {
                event_sum += event->prob;
            }
        }

        for (std::list<Event*>::iterator it = events_.begin();
             it != events_.end(); ++it) {
            (*it)->prob /= event_sum;
        }

        
        // Add prior contribution to ll
        float prior_ll = 0.0;
        
        for (std::list<Event*>::iterator it = events_.begin();
             it != events_.end(); ++it) {
            Event* event = *it;
                // Update
            prior_ll += log(event->prob + 1.11e-16);   
            
        }
        prior_ll = alpha_ * prior_ll;
        ll -= prior_ll;
                
        delta_loglik = abs(old_loglik - ll);
        old_loglik = ll;
    }    
}

void
NucPosFinder::UpdateEmpiricalDist()
{
    Eigen::VectorXd dist = Eigen::VectorXd::Zero(width_);
    
    for(std::list<Event*>::iterator it = events_.begin();
        it != events_.end(); ++it) {
        Event* event = *it;
        int start = (int)std::max(0, event->pos - width_);
        int stop = (int)std::min(pileup_->stop(), event->pos + width_);
        for (int i = start; i < stop; ++i) {
            int temp = pileup_->get(i);
            dist(abs(event->pos - i)) += temp;
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