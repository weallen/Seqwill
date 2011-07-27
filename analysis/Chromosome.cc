#include "analysis/Chromosome.h"


bool LoadChr(const std::string& fname, const std::string& genome_name,
            const std::string& chrname, Chromosome::Ptr chr)
{
  TrackFile trackio;
  trackio.Open(fname);
  Track<unsigned char>::Ptr track(new Track<unsigned char>);

  if (trackio.ReadSubTrack<unsigned char>(genome_name, chrname, track)) {
    chr->set_data(*track);
    chr->set_name(chrname);
    return true;
  }
  return false;
}

//
// CHROMOSOME
//
bool SaveChrFromFASTA(const std::string& outname, const std::string& seqname,
                      const std::string& genome_name)
{
  // parse file
  std::string currseq;
  FlatFileParser fastaparse;
  std::string line;
  TrackFile trackio;
  trackio.Open(outname);
  Track<unsigned char>::Ptr track(new Track<unsigned char>);

  DEBUGLOG("Saving chr " + seqname + " for " + genome_name);


  fastaparse.Open(seqname);

  // count number of lines in fasta
  std::string chrname;
  while (fastaparse.GetLine(line)) {
    // skip the header
    if (line[0] == '>')
      chrname = line.substr(1);
    else {
      currseq.append(line);
    }
  }
  track->set_subtrackname(chrname);
  track->set_extends(0, currseq.length());
  std::copy(currseq.begin(), currseq.end(), track->begin());
  bool ret = trackio.WriteSubTrack<unsigned char>(genome_name, track);
  trackio.Close();
  return ret;
}
