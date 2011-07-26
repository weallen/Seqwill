#include "analysis/Genome.h"


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

Track<float>::Ptr
GenomeData::GetTrackForChrom(const std::string& trackname, const std::string& chrname)
{
  Track<float>::Ptr curr_track;
  ChrMap m = open_chrs_[trackname];
  if (m.find(chrname) != m.end()) {
    Track<float>::Ptr curr_track = m[chrname];    
  } else {
    Track<float>::Ptr curr_track(new Track<float>);    
    trackfile_->ReadSubTrack<float>(trackfile_name_, trackname, chrname, curr_track);
    m[chrname] = curr_track;
  }
  return curr_track;
}

void
GenomeData::SaveTrackFromWIG(const std::string& wigname, const std::string& trackname, int resolution)
{
  WIGFile wig;
  
  DEBUGLOG("Loading track " + trackname + " from WIG file " + wigname);    
  ParseWig(wigname, &wig);
  std::vector<std::string> chrnames = genome_info_.chr_names();
  
  for (std::vector<std::string>::const_iterator it = chrnames.begin();
       it != chrnames.end(); ++it) {
    Track<float>::Ptr track(new Track<float>);
    track->set_resolution(resolution);
    track->set_extends(0, genome_info_.chr_size(*it));
    track->set_trackname(trackname);
    track->set_subtrackname(*it);
    
    std::vector<WIGLine> curr_chr = wig.lines[ChrToNum(*it)];
    for (std::vector<WIGLine>::const_iterator it = curr_chr.begin();
         it != curr_chr.end(); ++it) {
      
    }
  }
}
