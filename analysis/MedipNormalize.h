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

struct Frag
{  
  int start;
  int stop;
  int* cpgs;
  float* cpg_probs;
  int num_cpgs;
};

struct CpG
{
  int pos;
  float weight;
};

class MedipNormalize : public Analysis<float, float>
{
public:
  using Analysis<float, float>::analysis_name_;
  
  MedipNormalize()    
  : frag_len_(-1) 
  , bin_size_(10000)
  { analysis_name_ = std::string("MedipNormalize"); }

  virtual ~MedipNormalize();

  void set_frag_len(int frag_len)
  { frag_len_ = frag_len; }
  
  void set_chr(Track<unsigned char>::Ptr chr)
  { chr_ = chr; }
  
  void set_bam(BamIO* b) 
  { bio_ = b; }
  
private:  
  void ReadsToFrags();  
  void FindCpG();
  void AssignCpGToFrags();
  
  int PosToBinIdx(int pos);
  
  virtual void ComputeAnalysis();
                         
  Track<unsigned char>::Ptr chr_;
  int frag_len_;
  std::vector<Frag> frags_;
  // size of bin to order CpGs by
  int bin_size_;
  std::vector<std::vector<CpG> > cpgs_;
  BamIO* bio_;
};



// divides genome into intervals and counts CpGs 
// (possibly and CpAs) in each interval
class CpGCounter : public Analysis<unsigned char, int>
{
public:
  CpGCounter() 
  : do_cpa_(false)
    , res_(1)
    , tname_("")
    , stname_("")
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
  std::string tname_;
  std::string stname_;
};
#endif
