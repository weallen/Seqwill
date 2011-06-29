#include "analysis/Genome.h"




// GenomeMgr

void Genome::GetChrSeq(const std::string& chrname, ChromosomePtr& chr)
{
  if (chr_file_ == "") {
    return -1;
  } else {

    LoadChr(fname, chr);
    return 1;
  }
}
