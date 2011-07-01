#include "analysis/Chromosome.h"

int SaveChrFromFASTA(const std::string& fname, const std::string& seqfname)
{
}

int LoadChr(const std::string& fname, const std::string& chrname, Chromosome* chr)
{
}

//
// CHROMOSOME
//
void SaveChrSeqFromFASTA(const std::string& fname, const std::string& seqfname)
{
  // parse file
  std::string currseq;
  FlatFileParser fastaparse;
  std::string line;
  TrackWriter<const char> writer;

  writer.Open(fname);
  fastaparse.Open(seqfname);
  // count number of lines in fasta
  int i = 0;
  std::string chrname = "";
  while (fastaparse.GetLine(line)) {
    // skip the header
    if (line[0] == '>')
      chrname = line.substr(1);
    else {
      currseq.append(line);
    }
  }
  writer.WriteSubTrack(chrname, (const char*) currseq.c_str());
  writer.Close();
}
