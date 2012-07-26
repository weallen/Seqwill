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
    assert(stop_ != 0);
    output_ = Track<int>::Ptr(new Track<int>);
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
            int pos = (int)floor((float)(it->second.partner_pos + it->second.pos + it->second.len)/2.0);
            float temp = output_->get(pos);
            output_->set(pos, temp+1.0);
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
 
void
ExtendPileup::ComputeProcess()
{
    assert(stop_ != 0);
    output_ = Track<int>::Ptr(new Track<int>);
    output_->set_trackname(tname_);
    output_->set_subtrackname(stname_);
    output_->set_extends(start_, stop_);

    for (size_t i = 0; i < output_->size(); ++i) {
        output_->set(i, 0);
    }

    for (SingleReadFactory::iterator it = reads_->begin();
         it != reads_->end(); ++it) {
      int start = -1;
      int stop = -1;

      if (it->second.strand == kFwd) {
        start = it->second.pos;
        stop = it->second.pos + fraglen_ - 1;
      } else if (it->second.strand == kRev) {
        start = it->second.pos + it->second.len - fraglen_;
        stop = it->second.pos + it->second.len - 1;
      }
      int pos = floor(((float)start + stop)/2.0);
      int temp = output_->get(pos);
      output_->set(pos, temp+1);
    }
}


//-------------------------------------------------

void
NucKDE::ComputeAnalysis()
{
    output_ = Track<float>::Ptr(new Track<float>);
    output_->set_extends(0, input_->size());
    output_->set_trackname(tname_);
    output_->set_subtrackname(stname_);

    Track<float>::Ptr temp(new Track<float>);

    temp->set_extends(0, input_->size());
    for (size_t i = 0; i < temp->size(); ++i) {
      temp->set(i, 0.0);
    }

    float* dists = new float[2*w_];
    
    int j = 0;
    for (int i = -w_; i < w_; ++i) {
      dists[j] = pow(1.0 - (i / (float)w_)*(i / (float)w_), 3);
      j++;
    }    

    tbb::parallel_for(tbb::blocked_range<size_t>(w_, (input_->size()-w_), 1000),
		      TBB_NucKernel(*input_, *temp, dists, w_));
    
   // std::cout << "Got here" << std::endl;
    
    tbb::parallel_for(tbb::blocked_range<size_t>(150, input_->size() - 150, 1000),
                      TBB_NucPositioning(*temp, *output_, w_));
    temp.reset();    
}

//-------------------------------------------------

void
TBB_NucKernel::operator()(const tbb::blocked_range<size_t>& r) const
{      
    for (size_t pos = r.begin(); pos != r.end(); ++pos) {        
        float val = D(pos);        
        output_.set(pos, val); 
    }
}



float 
TBB_NucKernel::D(int i) const
{
    float val = 0.0;
    for (int j = 0; j < 2 * w_; ++j) {
        val += dists_[j] * (float)track_.get(i - w_ + j);
    }
    return val;
}

//-------------------------------------------------

void
TBB_NucPositioning::operator()(const tbb::blocked_range<size_t>& r) const
{
    for (size_t pos = r.begin(); pos != r.end(); ++pos) {        
        float val = temp_.get(pos);
        int start = pos - 150;
        int stop = pos + 150;
        float denom = 0.0;
        for (int i = start; i < stop; ++i) {
            denom +=  temp_.get(i);
        }
        denom *= (1.09 / w_);   
	if (denom > 0) {
	  output_.set(pos, val / denom);
	}
    }
}

//-------------------------------------------------




void
NucPositioner::FindMaxima()
{
  for (size_t i = 5; i < input_->size()-5; ++i) {
    if (input_->get(i) >= input_->get(i+1) 
	&& input_->get(i) >= input_->get(i-1)
	&& input_->get(i) >= input_->get(i+2)
	&& input_->get(i) >= input_->get(i-2)
	&& input_->get(i) >= input_->get(i+3)
	&& input_->get(i) >= input_->get(i-3) &&
	input_->get(i) > 0.0) {
      Nuc n;
      n.pos = i;
      n.weight = input_->get(i);
      nuc_maxima_.push_back(n);
    }
  }
}

void
NucPositioner::SelectNucs()
{
  nucs_.clear();
  int width = floor(nuc_size_/2.0);
  size_t i = 0;
  while (i < nuc_maxima_.size()-1) {
    if (nuc_maxima_[i].pos + width + 10 > nuc_maxima_[i+1].pos - width) {
      if (nuc_maxima_[i].weight > nuc_maxima_[i+1].weight) {
	nucs_.push_back(nuc_maxima_[i]);
      } else {
	nucs_.push_back(nuc_maxima_[i+1]);
      }
      i += 2;
    } else {
      nucs_.push_back(nuc_maxima_[i]);
      i++;
    }
  }
}

std::vector<BEDelement>
NucPositioner::NucPosAsBED()
{
  std::vector<BEDelement> elems;
  int width = floor(nuc_size_/2.0);
  for (std::vector<Nuc>::iterator it = nucs_.begin();
       it != nucs_.end(); ++it) {
    elems.push_back(BEDelement(chrname_, 
			       it->pos - width, it->pos + width, 
			       Stringify(it->weight)));    
  }
  return elems;
}
