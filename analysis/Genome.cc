#include "analysis/Genome.h"




// GenomeMgr

int Genome::GetChrSeq(const std::string& chrname, Chromosome* chr)
{
  if (chr_file_ == "") {
    return -1;
  } else {

    LoadChr(chr_file_, chrname, chr);
    return 1;
  }
}
