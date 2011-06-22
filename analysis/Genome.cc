#include "analysis/Genome.h"

//
// GENOME
//

void Genome::Open(const char* dirname) 
{
  std::string s(fname);
  Open(s);
}

void Genome::Open(const std::string& dirname) 
{ 
  m_dirname = dirname;
  if (!m_isopen) {
    m_data_dir = opendir(m_dirname.c_str());
    if (m_data_dir == NULL) {
      mkdir(m_dirname.c_str(), S_IRWXU); 
      m_data_dir = opendir(m_dirname.c_str());
    } else {
      struct dirent* dirconts; 
      while((dirconts = readdir(m_data_dir))) {
        string s(dirconts->d_name);
        string extension(".h5");
        if (Contains(s, extension)) {
          // take all but last 3 chars: e.g. "chr1.h5" -> "chr1"
          m_chrnames.push_back(s.substr(0, s.size() - 3));
        }
      }
    }
    m_isopen = true;
  }
}

void Genome::Close() 
{
  if (m_isopen) {
    for (map<string, Chromosome*>::iterator i = m_open_chrs.begin();
      i != m_open_chrs.end(); ++i) {
      (*i).second->Close();
    }

    if (m_data_dir != NULL)
      closedir(m_data_dir);
    m_isopen = false;
  }
}


void Genome::LoadSeq(const string& seqname) 
{
  // parse file
  string currseq;
  FlatFileParser fastaparse;
  string line;
  fastaparse.Open(seqname);
  // count number of lines in fasta
  int i = 0;
  string chrname = "";
  while (fastaparse.GetLine(line)) {
    // skip the header
    if (line[0] == '>') 
      chrname = line.substr(1);
    else {
      currseq.append(line);
    }
  }
  Chromosome chr(chrname, m_dirname);
  chr.WriteSeq(currseq);
  m_chrnames.push_back(chrname);
}

// Chromosome stuff
svec<string> Genome::GetAllChromosomeNames() 
{
  if (!m_isopen) {
    Open();
  }
  return m_chrnames;
}

bool Genome::HasChromosome(const string& chrname) {
  for (svec<string>::iterator i = m_chrnames.begin(); 
      i != m_chrnames.end(); ++i) {
    if (chrname == *i) 
      return true;
  }
  return false;
}

Chromosome* Genome::GetChromosome(const std::string& chrname) 
{
  Chromosome* chrout;
  if (!HasChromosome(chrname)) {
    return -1;
  }
  if (m_open_chrs.find(chrname) == m_open_chrs.end()) {
    Chromosome* chr = new Chromosome(chrname, m_dirname);
    chr->Open();
    m_open_chrs[chrname] = chr;
  }
  chrout = m_open_chrs[chrname];
  return chrout;
}

