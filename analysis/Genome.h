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

#define SUPERCONTIG_LEN 100000

using namespace boost::numeric;

class Genome
{
 public:
  Genome()
    : m_dirname("")
    , m_data_dir(NULL)
    , m_isopen(false) {}
  
  virtual ~Genome() { 
    Close();
    for (map<string, Chromosome*>::iterator i = m_open_chrs.begin();
	 i != m_open_chrs.end(); ++i) {
      delete (*i).second;
    }
  }

    void Open(const char* dirname);   
    void Open(const std::string& dirname);
    void Close();
    bool IsOpen() { return m_isopen; }

    // Sequence stuff
    void LoadSeq(const std::string& seqname);

    // Chromosome stuff
    std::vector<string> GetAllChromosomeNames();
    bool HasChromosome(const std::string& chrname);
    Chromosome* GetChromosome(const std::string& chrname);  

    std::vector<Chromosome*> GetChromosomes() { return m_chrs; }
    // Track stuff
    svec<string> GetAllTrackNames() {
      if (!m_isopen) {
        Open();
      }            
    }

  private: 
    void LoadTrackNames();

    string m_dirname;
    DIR* m_data_dir;
    bool m_isopen;

    svec<string> m_chrnames;
    std::map<std::string, Chromosome*> m_open_chrs;
    std::vector<Chromosome*> m_chrs;
    svec<string> m_tracknames;
};

#endif
