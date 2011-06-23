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

#include <boost/numeric/ublas/vector.hpp>

#include <hdf5.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "analysis/Track.h"
#include "analysis/Chromosome.h"
#include "hdf5/HDFFile.h"

#define SUPERCONTIG_LEN 100000

class Genome
{
 public:
  Genome()
    : dirname_("")
    , data_dir_(NULL)
    , isopen_(false) {}
  
  virtual ~Genome() { 
    Close();
    for (map<string, Chromosome*>::iterator i = open_chrs_.begin();
	 i != open_chrs_.end(); ++i) {
      delete (*i).second;
    }
  }

    void Open(const char* dirname);   
    void Open(const std::string& dirname);
    void Close();
    bool IsOpen() { return isopen_; }

    // Sequence stuff
    void LoadSeq(const std::string& seqname);

    // Chromosome stuff
    std::vector<string> GetAllChromosomeNames();
    bool HasChromosome(const std::string& chrname);
    Chromosome* GetChromosome(const std::string& chrname);  

    std::vector<Chromosome*> GetChromosomes() { return chrs_; }
    
    // Track stuff
    svec<string> GetAllTrackNames() {
      if (!isopen_) {
        Open(dirname_);
      }
      return tracknames_;
    }

  private: 
    void LoadTrackNames();

    string dirname_;
    DIR* data_dir_;
    bool isopen_;

    svec<string> chrnames_;
    std::map<std::string, Chromosome*> open_chrs_;
    std::vector<Chromosome*> chrs_;
    svec<string> tracknames_;
};

#endif
