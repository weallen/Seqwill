#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <string>

#include <boost/numeric/ublas/vector.hpp>

//#include <H5Cpp.h>
#include <hdf5.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "analysis/Track.h"

using namespace std;
using namespace boost::numeric;

class Chromosome
{
  public:
    Chromosome(const string& chrname, const string& dirname)
      : m_dirname(dirname)
      , m_isopen(false)
      , m_chrname(chrname)
      , m_len(-1)
    { }

    virtual ~Chromosome() { }

    void Open();
    int Close();

    bool IsOpen() { return m_isopen; }

    const string& GetName() { return m_chrname; }
    int WriteSeq(const string& seq);
    int ReadSeq(string* seq);

    int GetLength();

    svec<string> GetTrackNames();
    int WriteTrack(const string& name, const ublas::vector<float>& v);
    ublas::vector<float> ReadTrack(const string& name);
    void DeleteTrack(const string& name);
    
  private:
    string m_dirname;
    bool m_isopen;
    string m_chrname;
    int m_len;
    hid_t m_h5file;
    svec<string> m_tracknames;
};

#endif
