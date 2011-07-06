#include "analysis/Chromosome.h"


int LoadChr(const std::string& fname, const std::string& chrname, Chromosome* chr)
{
  return 1;
}

//
// CHROMOSOME
//
void SaveChrFromFASTA(const std::string& outname, const std::string& seqname,
                      const std::string& genome_name)
{
  // parse file
  std::string currseq;
  FlatFileParser fastaparse;
  std::string line;
  TrackIO::Ptr trackio(new TrackIO);
  typename Track<char>::Ptr track(new Track<char>);

  DEBUGLOG("Saving chr " + seqname + " for " + genome_name);


  trackio->Open(outname);
  fastaparse.Open(seqname);

  // count number of lines in fasta
  int i = 0;
  std::string chrname;
  while (fastaparse.GetLine(line)) {
    // skip the header
    if (line[0] == '>')
      chrname = line.substr(1);
    else {
      currseq.append(line);
    }
  }
  track->set_name(chrname);
  track->set_extends(0, currseq.length());
  std::copy(currseq.begin(), currseq.end(), track->begin());
  trackio->WriteSubTrack<char>(chrname, track);
}
