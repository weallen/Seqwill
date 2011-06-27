#ifndef SEQDATA_H_
#define SEQDATA_H_
#include <string>
#include "base/DNASequence.h"

class SeqData 
{
public:
  std::string chr_name;
  int start;
  int end;
  DNASequence seq;
};
#endif
