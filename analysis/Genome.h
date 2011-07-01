#ifndef GENOME_H_
#define GENOME_H_

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include <hdf5.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"

#include "io/TrackWriter.h"
#include "io/TrackReader.h"

#include "analysis/Chromosome.h"

class Genome
{
public:
  Genome() {}
  virtual ~Genome() {}

  void SetChrFile(const std::string& fname) { chr_file_ = fname; }


  int GetChrSeq(const std::string& chrname, Chromosome* chr);

private:

  std::map<std::string, Chromosome*> open_chrs_;
  std::string chr_file_;
};


#endif
