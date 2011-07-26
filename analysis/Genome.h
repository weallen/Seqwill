#ifndef GENOME_H_
#define GENOME_H_

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include <boost/shared_ptr.hpp>

#include "base/Types.h"
#include "base/Log.h"
#include "base/Common.h"
#include "base/WIG.h"
#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/Types.h"
#include "io/TrackIO.h"

class GenomeInfo;
class GenomeData;

void SaveTrackFromWIG(const std::string& wigname, const std::string& trackname);
void LoadGenomeInfoFromChr(const std::string& chrtracksname, const std::string& genomename, GenomeInfo* info);

class GenomeInfo 
{
public:
  GenomeInfo() {}
  virtual ~GenomeInfo() { }
    
  const std::vector<std::string> chr_names() const
  { return chr_names_; }
  
  const int chr_size(const std::string& chrname) 
  { return chr_sizes_[chrname]; }
  
  void set_chr_names(const std::vector<std::string>& chrnames)
  { chr_names_ = chrnames; }
  
  void set_chr_size(const std::string& chrname, int size)
  { chr_sizes_[chrname] = size; }
  
private:
  DISALLOW_COPY_AND_ASSIGN(GenomeInfo)
  
  std::vector<std::string> chr_names_;
  std::map<std::string, int> chr_sizes_;
};

// GenomeData interfaces between TrackFile and user.
// Holds a collection of tracks -- one for each chromosome
// And information about the genome itself
class GenomeData
{
public:
  typedef std::map<std::string, Track<float>::Ptr> ChrMap;  
  typedef std::map<std::string, ChrMap> SubtrackToChrMap;
  
  GenomeData() {}
  virtual ~GenomeData() {}
  
  void set_genome_info(const GenomeInfo& ginfo);
  
  const GenomeInfo& genome_info() const 
  { return genome_info_; }
                     
  Track<float>::Ptr GetTrackForChrom(const std::string& trackname, const std::string& chrname);
  
  void SaveTrackFromWIG(const std::string& wigname, const std::string& trackname, int resolution);
  
  void InitEmpty();
  
private:
  std::string trackfile_name_;
  TrackFile::Ptr trackfile_;  
  GenomeInfo genome_info_;
  SubtrackToChrMap open_chrs_;
};


#endif
