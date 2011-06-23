#include "analysis/Genome.h"

//
// GENOME
//

void Genome::Open(const char* dirname) 
{
  std::string s(dirname);
  Open(s);
}

void Genome::Open(const std::string& dirname) 
{ 
  dirname_ = dirname;
  if (!isopen_) {
    data_dir_ = opendir(dirname_.c_str());
    if (data_dir_ == NULL) {
      mkdir(dirname_.c_str(), S_IRWXU); 
      data_dir_ = opendir(dirname_.c_str());
    } else {
      struct dirent* dirconts; 
      while((dirconts = readdir(data_dir_))) {
        string s(dirconts->d_name);
        string extension(".h5");
        if (Contains(s, extension)) {
          // take all but last 3 chars: e.g. "chr1.h5" -> "chr1"
          chrnames_.push_back(s.substr(0, s.size() - 3));
        }
      }
    }
    isopen_ = true;
  }
}

void Genome::Close() 
{
  if (isopen_) {
    for (std::map<std::string, Chromosome*>::iterator i = open_chrs_.begin();
      i != open_chrs_.end(); ++i) {
      (*i).second->Close();
    }

    if (data_dir_ != NULL)
      closedir(data_dir_);
    isopen_ = false;
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
  Chromosome chr(chrname, dirname_);
  chr.WriteSeq(currseq);
  chrnames_.push_back(chrname);
}

// Chromosome stuff
svec<string> Genome::GetAllChromosomeNames() 
{
  if (!isopen_) {
    Open();
  }
  return chrnames_;
}

bool Genome::HasChromosome(const string& chrname) {
  for (svec<string>::iterator i = chrnames_.begin(); 
      i != chrnames_.end(); ++i) {
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
  if (open_chrs_.find(chrname) == open_chrs_.end()) {
    Chromosome* chr = new Chromosome(chrname, dirname_);
    chr->Open();
    open_chrs_[chrname] = chr;
  }
  chrout = open_chrs_[chrname];
  return chrout;
}

