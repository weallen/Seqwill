#include "analysis/MedipNormalize.h"



/* ------------------------------ */
MedipNormalize::~MedipNormalize() {
  if (cpgs_ != NULL)
    delete[] cpgs_;
}

void
MedipNormalize::ReadsToFrags()
{
  assert(reads_ != NULL);

  int nfrags = 0;
  // index frags by position of first mate
  if (frag_len_ > 0) {
    for (SingleReadFactory::iterator it = reads_->begin();
         it != reads_->end(); ++it) {
      Frag f;
      f.num_cpgs = 0;
      if (it->second.strand == kFwd) {
        f.start = it->second.pos;
        f.stop = it->second.pos + frag_len_ - 1;
      } else if (it->second.strand == kRev) {
        f.start = it->second.pos + it->second.len - frag_len_;
        f.stop = it->second.pos + it->second.len - 1;
      }
      frags_.push_back(f);
      nfrags++;
    }    
  } else {
    for (SingleReadFactory::iterator it = reads_->begin();
         it != reads_->end(); ++it) {
      Frag f;
      f.num_cpgs = 0;
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
  }
  num_frags_ = nfrags;
}

void
MedipNormalize::FindCpG()
{  
  // initialize bins  
  std::vector<CpG> cpgs;
  
  int ncpgs = 0;
  for (size_t i = 0; i < chr_->size() - 1; ++i) {
    unsigned char curr = chr_->get(i);
    unsigned char next = chr_->get(i+1);
    if ((curr == 'C' && next == 'G')
        || (curr == 'c' && next == 'g')) {
      CpG cpg;
      cpg.pos = i;
      cpg.weight = 1.0;
      cpgs.push_back(cpg);
      ncpgs++;
    }
  }
  cpgs_ = new CpG[cpgs.size()];    
  std::copy(cpgs.begin(), cpgs.end(), &(cpgs_[0]));
  num_cpgs_ = ncpgs;
}

int
MedipNormalize::PosToBinIdx(int pos)
{
  return floor((float)pos / bin_size_);
}

void
MedipNormalize::AssignCpGToFrags()
{
  
  // index CpGs
  num_bins_ = ceil((float)chr_->size() / bin_size_);
  std::vector<std::vector<int> > idx(num_bins_+1);
  for (int i = 0; i < num_cpgs_; ++i) {
    idx[PosToBinIdx(cpgs_[i].pos)].push_back(i);
  }
  num_bins_ = ceil(chr_->asize() / (float)bin_size_);
  std::vector<int> index(num_bins_);
  //std::cerr << index.size() << std::endl;
  for (int frag_idx = 0; frag_idx < num_frags_; ++frag_idx) {   
    Frag f = frags_[frag_idx];
    f.num_cpgs = 0;
    std::vector<CpG*> frag_cpgs;
    int start = std::max(PosToBinIdx(f.start)-1,0);
    int stop = std::min(PosToBinIdx(f.stop)+1, num_bins_);
    
    for (int i = start; i <= stop; ++i) { 
      if (idx[i].size() > 0) {
        std::vector<int> cpg_idx = idx[i];
        for (size_t j = 0; j < cpg_idx.size(); ++j) {
          int curr_idx = cpg_idx[j];
          if ((cpgs_[curr_idx].pos >= f.start) && (cpgs_[curr_idx].pos <= f.stop)) {
            frag_cpgs.push_back(&(cpgs_[curr_idx]));
          }
        }
      } 
    }
    f.num_cpgs = (int)frag_cpgs.size();
    for (size_t k = 0; k < frag_cpgs.size(); ++k) {        
      f.cpg_probs.push_back(1.0 / (float)frag_cpgs.size());
      f.cpgs.push_back(frag_cpgs[k]);
    }
    frags_[frag_idx] = f;
  }    
}

void
MedipNormalize::IterativelyReweightCpGs()
{
  float delta_loglik = INFINITY;
  float old_loglik = 0.0;
  int n = 0;

  while (delta_loglik > 0.01) {
    n++;
    //DEBUGLOG("REWEIGHT STEP " + Stringify(n));
    // zero all the CpGs
    for (int i = 0; i < num_cpgs_; ++i) {
      cpgs_[i].weight = 0.0;
    }
    
    // Update weights given current per fragment cpg probs
    for (int i = 0; i < num_frags_; ++i) {
      for (size_t j = 0; j < frags_[i].cpgs.size(); ++j) {
        frags_[i].cpgs[j]->weight += frags_[i].cpg_probs[j];
      }
   }
 
    // Update per frag cpg probs by weights
    float ll = 0.0;
   
    for (int i = 0; i < num_frags_; ++i) {
      float frag_ll = 0.0;
      float sum = 0.0;
      for (int j = 0; j < frags_[i].num_cpgs; ++j) {
        float temp = frags_[i].cpgs[j]->weight;
        frags_[i].cpg_probs[j] = temp;
        sum += temp;
        ll += log(temp+1.11e-16);
      }
      // normalize to make prob dist over cpgs
      for (int j = 0; j < frags_[i].num_cpgs; ++j) {
        frags_[i].cpg_probs[j] /= sum;
      }
    }
    delta_loglik = abs(old_loglik - ll);
    old_loglik = ll;
    //    std::cerr << ll << std::endl;
  }
}

void
MedipNormalize::ComputeAnalysis()
{
//  assert(!isnan(beta_) && !isnan(alpha_) && frag_len_ != -1);
  //DEBUGLOG("Finding CpGs " + stname_ + "...");
  FindCpG();
  //DEBUGLOG("Assigning reads to frags " + stname_ + "...");
  ReadsToFrags();
  //DEBUGLOG("Assinging cpgs to frags " + stname_ + "...");
  AssignCpGToFrags();
  DEBUGLOG("Iteratively reweighting cpgs " + stname_ + "...");
  IterativelyReweightCpGs();
  
  output_ = TrackOutPtr(new Track<float>);
  output_->set_resolution(res_);
  output_->set_extends(0, floor(((float)chr_->stop())/res_));
  output_->set_trackname(tname_);
  output_->set_subtrackname(stname_);
  for (size_t t = 0; t < output_->size(); ++t) {
    output_->set(t, 0.0);
  }
  for (int i = 0; i < num_cpgs_; ++i) {
    int bin_idx = floor((float)cpgs_[i].pos / res_);
    float curr = output_->get(bin_idx);
    output_->set(bin_idx, curr + cpgs_[i].weight);
  }
}



//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void 
CpGCounter::ComputeAnalysis() 
{
  assert(tname_ != "" && stname_ != "");
  output_ = TrackOutPtr(new Track<int>);
  output_->set_resolution(res_);
  output_->set_abs_extends(0, input_->stop());
  output_->set_extends(0, floor(((float)input_->stop())/res_));
  output_->set_trackname(tname_);
  output_->set_subtrackname(stname_);
  for (size_t i = 0; i < output_->size(); ++i) {
    output_->set(i, 0);
  }
  
  
  for (size_t i = 0; i < input_->size() - 1; ++i) {
    int curr_bin = floor((float) i / res_);
    unsigned char curr = input_->get(i);
    unsigned char next = input_->get(i+1);
    if ((curr == 'C' && next == 'G')
        || (curr == 'c' && next == 'g')) {
      int old = output_->get(curr_bin);
      output_->set(curr_bin, old+1);
    }
  }
}
