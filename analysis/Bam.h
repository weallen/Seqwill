#ifndef BAM_H_
#define BAM_H_

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>

#include <bamtools/api/BamReader.h>

#include "analysis/Genome.h"


namespace boost::numeric::ublas = ublas;


class BamLoader
{
 public:

  BamLoader()  {}

  virtual ~BamLoader();

  // creates index, inits data structures
  void Init(Genome* g, const std::string& fname);

  void Close();
  
  void ReadBamFile();
  
  Genome* genome() { return genome_; }

  int extend() { return extend_; }
  
  void extend(int e) { extend_ = e; }
  
 private:
  void InitChrs();
  void PairedRead(const BamAlignment& read);
  void UnpairedRead(const BamAlignment& read);
  
  std::map<std::string, ublas::vector<float> > chrs_;
  BamReader reader_;  
  Genome* genome_;
  int extend_;
};


#endif
