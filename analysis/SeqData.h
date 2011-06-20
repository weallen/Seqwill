#ifndef DATA_H_
#define DATA_H_

#include <string>

#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamIndex.h>
#include <bamtools/api/BamAlignment.h>

#include "analysis/Genome.h"

using namespace std;

class SeqDataLoader
{
  public:
    SeqDataLoader() {}

    SeqDataLoader(const string& bamfilename) 
      : m_bamfilename(bamfilename)
      , m_reader()
    { }

    virtual ~SeqDataLoader() {}

  private:
    string m_bamfilename;
    BamTools::BamReader m_reader;
};

class SeqData
{
  public:
    SeqData() {}
    virtual ~SeqData() {}

      
  private:
    Genome m_genome;
};

#endif 
