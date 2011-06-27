#ifndef BAM_H_
#define BAM_H_

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

#include "analysis/Genome.h"
#include "base/StringUtil.h"
#include "hdf5/HDFFile.h"

enum Stand {
  FWD = 0,
  REV = 1,
  UNKOWN = 2
};

class ReadHit 
{
  uint64_t insert_id;
  uint64_t ref_id;
  int pos;

  uint64_t partner_ref;
  int partner_pos;

  Strand strand;
  // from cigar string
  int len;
};

class MateHit
{
 public:
  MateHit(const ReadHit* left_align,
	  const ReadHit* right_align)
    : left_align_(left_align)
    , right_align_(right_align)
  {}

 private:
  const ReadHit* left_align_;
  const ReadHit* right_align_;
};

// Map between hashed numeric ref_ids and string refids
class HitTable
{
 public:
  HitTable() {}
  virtual ~HitTable() {}
  

 private:  
  inline uint64_t HashString(const char* __s) 
  {    
    uint64_t hash = 0xcbf29ce484222325ull;
    for ( ; *__s; ++__s) {
      hash *= 1099511628211ull;
      hash ^= *__s;
    }
    return hash;
  }   

  // member vars

};

class BamLoader
{
 public:

  BamLoader()  {}
  
  virtual ~BamLoader();

  // creates index, inits data structures
  void Init(GenomeMgr* g, const std::string& fname);

  void Close();
  
  void ReadBamFile();
  
  GenomeMgr* genome() { return genome_; }
  
 private:
  std::map<std::string, int> GetRefSeqInfo();
  void InitChrs();
  void PairedRead(const BamTools::BamAlignment& read);
  void UnpairedRead(const BamTools::BamAlignment& read);

  HitTable table_;
  std::map<std::string, std::vector<BamRead> > chrs_;
  BamTools::BamReader reader_;  
  std::map<std::string, int> chrlens_;
  GenomeMgr* genome_;
};


#endif
