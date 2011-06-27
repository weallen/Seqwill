#ifndef GENOMESEQ_H_
#define GENOMESEQ_H_
#include <string>
#include "base/DNASequence.h"

class GenomeSeq
{
public:
  GenomeSeq() {}

  virtual ~GenomeSeq() {}

  int start;
  int end;
  std::string chr;
  DNASequence seq;
};


#endif
