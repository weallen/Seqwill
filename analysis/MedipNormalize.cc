#include "analysis/MedipNormalize.h"


MedipNormalize::~MedipNormalize() {
    delete bio_; 
    for (size_t i = 0; i < frags_.size(); ++i) {
      delete[] frags_[i].cpgs;
      delete[] frags_[i].cpg_probs;
    }
}

void
MedipNormalize::ReadsToFrags()
{
  SingleReadFactory* reads = bio_->LoadChrSingleReads(chr_->subtrackname());
  // index frags by position of first mate
  
  for (SingleReadFactory::iterator it = reads->begin();
       it != reads->end(); ++it) {
    if (it->second.is_first) {  
      Frag f;
      if (frag_len_ > 0) {
        if (it->second.strand == kFwd) {
          f.start = it->second.pos;
          f.stop = it->second.pos + frag_len_ - 1;
        } else if (it->second.strand == kRev) {
          f.start = it->second.pos + it->second.len - frag_len_;
          f.stop = it->second.pos + it->second.len - 1;
        }
      } else {
        //XXX may need to fix this
        if (it->second.partner_ref == it->second.ref_id) {
          if (it->second.pos > it->second.partner_pos) {
            f.start = it->second.pos;
            f.stop = it->second.partner_pos;
          } else {
            f.start = it->second.partner_pos;
            f.stop = it->second.pos;
          }          
        }
      }
      frags_.push_back(f);
    }
  }
}

void
MedipNormalize::FindCpG()
{  
  // initialize bins  
  cpgs_ = std::vector<std::vector<CpG> >(ceil((float)chr_->asize() / bin_size_));  
  
  for (size_t i = 0; i < input_->size() - 1; ++i) {
    unsigned char curr = chr_->get(i);
    unsigned char next = chr_->get(i+1);
    if ((curr == 'C' && next == 'G')
        || (curr == 'c' && next == 'g')) {
      CpG cpg;
      cpg.pos = i;
      cpg.weight = 1.0;
      cpgs_[PosToBinIdx(i)].push_back(cpg);
    }
  }
}

void
MedipNormalize::AssignCpGToFrags()
{
  for (std::vector<Frag>::iterator it = frags_.begin();
       it != frags_.end(); ++it) {
    
  }    
}


int
MedipNormalize::PosToBinIdx(int pos) 
{
  return floor(((float)pos) / bin_size_);
}

void
MedipNormalize::ComputeAnalysis()
{
//  assert(!isnan(beta_) && !isnan(alpha_) && frag_len_ != -1);
  DEBUGLOG("Finding CpGs...");
  FindCpG();
  DEBUGLOG("Fitting linear regression...");
//  FitLinear();  
}



//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void 
CpGCounter::ComputeAnalysis() 
{
  assert(tname_ != "" && stname_ != "");
  output_ = TrackOutPtr(new Track<int>);
  output_->set_resolution(res_);
  output_->set_extends(floor(((float)input_->start()) / res_), floor(((float)input_->stop())/res_));
  output_->set_trackname(tname_);
  output_->set_subtrackname(stname_);
  unsigned char curr_char; 
  unsigned char prev_char;
  for (size_t i = 0; i < output_->size(); ++i) {
    output_->set(i, 0);
  }
  
  int curr_count;
  for (size_t i = 0; i < input_->size(); i += res_) {
    curr_count = 0;
    prev_char = input_->get(i);
    for (int j = 1; j < res_; ++j) {
      curr_char = input_->get(i+j);
      if (curr_char == 'G' && prev_char == 'C') 
	curr_count++;
      if (do_cpa_) {
	if (curr_char == 'A' && prev_char == 'C') {
	  curr_count++;
	}
      }
      prev_char = curr_char;
    }
    output_->set(i, curr_count);
  }
}
