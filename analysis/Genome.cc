#include "analysis/Genome.h"


void
SaveTracksFromWIG(const std::string& wigname, const std::string& trackfilename)
{
  TrackFile::Ptr trackfile(new TrackFile());
  std::vector<WIGLine> lines;
  
  DEBUGLOG("Loading track " + trackfilename + " from WIG file " + wigname);  
  ParseWig(wigname, &lines);
  std::vector<WIGLine>::const_iterator it;
  for (it = lines.begin(); it != lines.end(); ++it) {
    
  }
}

void 
LoadGenomeInfoFromChr(const std::string& chrtracksname, const std::string& genomename, GenomeInfo* info)
{
  TrackFile f;
  TrackMetadata m;
  int stop = -1;
  int resolution = 1;
  f.Open(chrtracksname);
  std::vector<std::string> chr_names = f.GetSubTrackNames(genomename);
  info->set_chr_names(chr_names);
  for (std::vector<std::string>::iterator it = chr_names.begin();
       it != chr_names.end(); ++it) {
    m = f.GetSubTrackMetadata(genomename, *it);
    stop = m.GetIntMetadata("Stop", -1);
    resolution = m.GetIntMetadata("Resolution", 1);
    info->set_chr_size(*it, stop * resolution);
  }
  f.Close();
}

//---------------------------------------------------------------------------

