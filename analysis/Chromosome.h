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
#include "base/Common.h"

class Chromosome 
{
public:
  typedef boost::shared_ptr<Chromosome> Ptr;
  typedef Track<unsigned char>::Ptr CharTrackPtr;

  Chromosome() {
   Init();
  }

  Chromosome(const std::string& chrname)
    : name_(chrname)
  {
    Init();
  }

  virtual ~Chromosome() {
  }

  void Init() {
    data_ = DNASequence::Ptr(new DNASequence);
  }

  std::string name() const { return name_; }
  int length() { return data_->size(); }

  DNASequence::Ptr data() { return data_; }

  void set_data(const Track<unsigned char>& data)
  {
    data_->Allocate(data.size());
    data_->Assign(data);
  }

  void set_name(const std::string& name) { name_ = name; }

private:
  DISALLOW_COPY_AND_ASSIGN(Chromosome)

  std::string name_;
  DNASequence::Ptr data_;
};

bool SaveChrFromFASTA(const std::string& outname, const std::string& seqname,
                         const std::string& genome_name);

bool LoadChr(const std::string& fname, const std::string& genome_name,
            const std::string& chrname, Chromosome::Ptr chr);

#endif
