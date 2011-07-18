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


//---------------------------------------------------------------------------

void 
GenomeInfo::Open(const std::string& fname)
{
}

void 
GenomeInfo::Load() 
{  
}

void
GenomeInfo::Close()
{
}

