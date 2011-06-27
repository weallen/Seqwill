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
#include "analysis/Track.h"
#include "analysis/Chromosome.h"
#include "hdf5/HDFTrackReader.h"
#include "hdf5/HDFTrackWriter.h"

#define SUPERCONTIG_LEN 100000

typedef boost::shared_ptr<GenomeMgr> GenomeMgrPtr;

int LoadGenomeMgr(const std::string& fname, GenomeMgr g);

class GenomeMgr
{
public:
  GenomeMgr() {}

  virtual ~GenomeMgr() {
    Close();
  }

  void Open(const char* fname);
  void Open(const std::string& fname);
  void Close();
  bool IsOpen() { return isopen_; }
  void Flush();

  // Sequence stuff
  void LoadChrSeq(const std::string& seqname);

  // Track stuff
  std::vector<string> GetChrNames();
  bool HasChromosome(const std::string& chrname);
  TrackPtr GetChromosome(const std::string& chrname, const std::string& trackname);
  TrackPtr GetChrRegion(const std::string& chrname, int start, int end, const std::string& trackname);

  void LoadTrackData(const std::string& tdfile);

  // Track stuff
  std::vector<string> GetTrackNames();

private:
  void LoadTrackNames();

  string filename_;
  bool isopen_;

  HDFTrackReader track_reader_;
  std::vector<string> chrnames_;
  std::map<std::string, int> chrlens_;
  std::vector<string> tracknames_;
};

#endif
