#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <string>

#include <boost/intrusive_ptr.hpp>

//#include <H5Cpp.h>
#include <hdf5.h>

#include "base/DNASequence.h"
#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/Log.h"
#include "base/RefCount.h"
#include "io/TrackIO.h"
#include "common/Track.h"


class Chromosome : public RefBase
{
public:
  typedef boost::intrusive_ptr<Chromosome> Ptr;

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

void SaveChrFromFASTA(const std::string& outname, const std::string& seqname,
                         const std::string& genome_name);
int LoadChr(const std::string& fname, const std::string& chrname, Chromosome::Ptr chr);

#endif
