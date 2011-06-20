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

using namespace std;
using namespace boost::numeric;

class Genome
{
  public:
    Genome() {}

    Genome(const char* dirname) 
      : m_dirname(dirname)
      , m_data_dir(NULL)
      , m_isopen(false)
    {}

    Genome(const string& dirname) 
      : m_dirname(dirname) 
      , m_data_dir(NULL)
      , m_isopen(false)
    {}


    virtual ~Genome() { 
      Close();
      for (map<string, Chromosome*>::iterator i = m_open_chrs.begin();
        i != m_open_chrs.end(); ++i) {
        delete (*i).second;
      }
    }

    void Open();   
    void Close();
    bool IsOpen() { return m_isopen; }

    // Sequence stuff

    void LoadSeq(const string& seqname);

    // Chromosome stuff
    svec<string> GetAllChromosomeNames();
    bool HasChromosome(const string& chrname);
    int GetChromosome(const string& chrname, Chromosome* chrout);  

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
    map<string, Chromosome*> m_open_chrs;

    svec<string> m_tracknames;
};

#endif
