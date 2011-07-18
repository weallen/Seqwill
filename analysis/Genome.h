#ifndef GENOME_H_
#define GENOME_H_

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sqlite3.h>

#include <boost/shared_ptr.hpp>

#include "base/Log.h"
#include "base/Common.h"
#include "base/WIG.h"
#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/Types.h"
#include "io/TrackIO.h"

class GenomeInfo;
class GenomeData;

void SaveGenomeDataFromWIG(const std::string& wigname, const std::string& genomefilename);

void LoadGenomeInfo(const std::string& infoname, GenomeInfo* g);


class GenomeInfo 
{
public:
  GenomeInfo() {}
  virtual ~GenomeInfo() { Close(); }
  
  void Open(const std::string& fname);
  void Close();
  
  const std::vector<std::string> chr_names() const
  { return chr_names_; }
  
  const int chr_size(const std::string& name) 
  { return chr_sizes_[name]; }
  
private:
  DISALLOW_COPY_AND_ASSIGN(GenomeInfo)
  void Load();
  void Save();
  
  std::string fname_;
  bool isopen_;
  std::vector<std::string> chr_names_;
  std::map<std::string, int> chr_sizes_;
};

class GenomeData
{
public:
  GenomeData() {}
  virtual ~GenomeData() {}
  
private:
  TrackFile::Ptr trackfile_;  
};


#endif
