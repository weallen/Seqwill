#ifndef MEDIP_NORMALIZE_H_
#define MEDIP_NORMALIZE_H_
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <math.h>

#include <Eigen/Dense>
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/SamSequenceDictionary.h>
#include <api/SamSequence.h>

#include "base/MatrixUtil.h"

#include "base/StringUtil.h"
#include "analysis/AnalysisBase.h"
#include "common/Track.h"
#include "analysis/Dist.h"
#include "io/BamIO.h"

/*
 * 1. set all cpgs to 1
 * repeat:
 *   2. assign fractions of fragments to CpGs based on scores
 *   3. sum fractions of reads for each CpG
 * for each fragment
 */


struct CpG
{
  int pos;
  float weight;
};

struct Frag
{  
  int start;
  int stop;
  std::vector<CpG*> cpgs;
  std::vector<float> cpg_probs;
  int num_cpgs;
};


class MedipNormalize : public Analysis<float, float>
{
public:
  using Analysis<float, float>::analysis_name_;
  
  MedipNormalize()    
  : frag_len_(-1) 
  , num_frags_(0)
  , bin_size_(10000)
  , num_cpgs_(0)
  , num_bins_(0)
  , cpgs_(NULL)
  , bio_(NULL)
  , res_(100)
  { analysis_name_ = std::string("MedipNormalize"); }

  virtual ~MedipNormalize();

  void set_frag_len(int frag_len)
  { frag_len_ = frag_len; }
  
  void set_chr(Track<unsigned char>::Ptr chr)
  { chr_ = chr; }
  
  void set_bam(BamIO* b) 
  { bio_ = b; }
  
  void set_bin_size(int bin) 
  { bin_size_ = bin; }
    
  const int num_cpgs() const
  { return num_cpgs_; }
  
  const int num_frags() const
  { return num_frags_; }
  
  std::vector<Frag>& frags() 
  { return frags_; }
  
  CpG* cpgs() 
  { return cpgs_; }
  
  void set_resolution(int res) 
  { assert(res > 0); res_ = res; }

  void ReadsToFrags();  
  void FindCpG();
  void AssignCpGToFrags();
  void IterativelyReweightCpGs();
  
protected:    
  virtual void ComputeAnalysis();
                        
  int PosToBinIdx(int pos);
  
  Track<unsigned char>::Ptr chr_;
  int frag_len_;
  std::vector<Frag> frags_;
  int num_frags_;
  // size of bin to order CpGs by
  int bin_size_;
  int num_cpgs_;
  int num_bins_;
  CpG* cpgs_;
  BamIO* bio_;
  int res_;
};



// divides genome into intervals and counts CpGs 
// (possibly and CpAs) in each interval
class CpGCounter : public Analysis<unsigned char, int>
{
public:
  CpGCounter() 
  : do_cpa_(false)
    , res_(1)
  { analysis_name_ = std::string("CpGCounter"); }
  
  virtual ~CpGCounter() {}
  
  void set_resolution(int res) 
  { assert(res > 0); res_ = res; }
  
  void set_do_cpa(bool t) 
  { do_cpa_ = t; }
 
  
private:
  virtual void ComputeAnalysis();
  bool do_cpa_;
  int res_;
};
#endif
