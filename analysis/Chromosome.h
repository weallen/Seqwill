#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <string>

#include <boost/shared_ptr.hpp>

//#include <H5Cpp.h>
#include <hdf5.h>

#include "base/DNASequence.h"
#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "io/TrackWriter.h"


class Chromosome;

int SaveChrFromFASTA(const std::string& fname, const std::string& seqfname);
int LoadChr(const std::string& fname, const std::string& chrname, Chromosome* chr);

class Chromosome
{
public:
  Chromosome() {}

  Chromosome(const std::string& chrname)
    : name_(chrname)
  { }

  virtual ~Chromosome() { }

  const std::string& GetName() { return name_; }

  int GetLength() { return len_; }

private:
  std::string name_;
  int len_;
  DNASequencePtr seq_;
};

#endif
