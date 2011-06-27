#include "analysis/Genome.h"

int LoadGenomeMgr(const std::string& fname, GenomeMgr g)
{
}

//
// GENOMEMGR
//

void GenomeMgr::Open(const char* dirname)
{
  std::string s(dirname);
  Open(s);
}

void GenomeMgr::Open(const std::string& dirname)
{ 
}

void GenomeMgr::Close()
{
  Flush();
}

// Deallocate everything
void GenomeMgr::Flush()
{
}

void GenomeMgr::LoadChrSeq(const string& seqname)
{
  track_reader_.Close();
  // parse file
  string currseq;
  FlatFileParser fastaparse;
  string line;
  HDFTrackWriter writer;

  writer.Open(filename_);

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
  chrnames_.push_back(chrname);
  writer.WriteChrSeq(chrname, currseq);
  writer.Close();
}

// Chromosome stuff
std::vector<string> GenomeMgr::GetChrNames()
{
  if (!isopen_) {
    Open();
  }
  return chrnames_;
}

bool GenomeMgr::HasChromosome(const string& chrname) {
  for (svec<string>::iterator i = chrnames_.begin(); 
      i != chrnames_.end(); ++i) {
    if (chrname == *i) 
      return true;
  }
  return false;
}

TrackPtr GenomeMgr::GetChromosome(const std::string& chrname, const std::string& trackname)
{
  TrackPtr track(new Track(trackname, chrname, 0, chrlens_[chrname]));
  TrackData td;
  track_reader_.GetChromosome();
}

TrackPtr GenomeMgr::GetChrRegion(const std::string& chrname, int start, int end, const std::string& trackname)
{

}
