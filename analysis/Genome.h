#ifndef GENOME_H_
#define GENOME_H_

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include <boost/shared_ptr.hpp>

#include "math.h"
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
  std::vector<std::string> chr_names_;
  std::map<std::string, int> chr_sizes_;
};

// GenomeData interfaces between TrackFile and user.
// Holds a collection of tracks -- one for each chromosome
// And information about the genome itself
class GenomeData
{
public:
  typedef std::map<std::string, Track<float>::Ptr > ChrMap;  
  
  GenomeData()
  : trackfile_name_("")
  , trackname_("")
  , trackfile_(new TrackFile()) 
  {}
  
  virtual ~GenomeData() { Close(); }
  
  void Init(const std::string& tfname_, const GenomeInfo& g); 
    
  void set_genome_info(const GenomeInfo& ginfo);
  void set_track_name(const std::string& trackname)
  { trackname_ = trackname; }
  
  const GenomeInfo& genome_info() const 
  { return genome_info_; }
                     
  Track<float>::Ptr GetSubTrackForChrom(const std::string& chrname);
  
  void SaveTrackFromWIG(const std::string& wigname, int resolution);
  
  void InitEmpty();
  
private:
  DISALLOW_COPY_AND_ASSIGN(GenomeData)
  
  void Close();
  
  std::string trackfile_name_;
  std::string trackname_;
  TrackFile::Ptr trackfile_;  
  GenomeInfo genome_info_;
  ChrMap open_chrs_;
};


#endif
