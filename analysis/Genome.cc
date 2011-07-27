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

void
GenomeData::Init(const std::string& fname, const GenomeInfo& g) 
{
  trackfile_name_ = fname;
  trackfile_->Open(fname);
  genome_info_ = g;
}

void
GenomeData::Close()
{
  trackfile_->Close();
}

 
Track<float>::Ptr
GenomeData::GetTrackForChrom(const std::string& trackname, const std::string& chrname)
{
  Track<float>::Ptr curr_track;
  ChrMap m = open_chrs_[trackname];
  if (m.find(chrname) != m.end()) {
    Track<float>::Ptr curr_track = m[chrname];    
  } else {
    Track<float>::Ptr curr_track(new Track<float>);    
    trackfile_->ReadSubTrack<float>(trackname, chrname, curr_track);
    m[chrname] = curr_track;
  }
  return curr_track;
}

void
GenomeData::SaveTrackFromWIG(const std::string& wigname, const std::string& trackname, int resolution)
{
  int curr_pos;
  WIGParser p;
  DEBUGLOG("Loading track " + trackname + " from WIG file " + wigname);    
  p.Open(wigname);
  std::vector<std::string> chrnames = genome_info_.chr_names();
  for (std::vector<std::string>::const_iterator it = chrnames.begin();
       it != chrnames.end(); ++it) {   
    std::cerr << *it << std::endl;
    Track<float>::Ptr track(new Track<float>);
    track->set_resolution(resolution);
    track->set_abs_extends(0, genome_info_.chr_size(*it));
    track->set_trackname(trackname);
    track->set_subtrackname(*it);
    track->set_extends(0, ceil(static_cast<float>(genome_info_.chr_size(*it))/resolution));
    bool ret = trackfile_->WriteSubTrack<float>(trackname, track);
    
  }
  
  //for (std::vector<WIGLine>::const_iterator j = curr_chr.begin();
//       j != curr_chr.end(); ++j) {
  //   curr_pos = j->pos;
//        track->aset(curr_pos, j->val);
  //    }
//      trackfile_->WriteSubTrack<float>(trackname, track);
    //} else {
     // ERRORLOG("Only FixedStep supported for now...");
    //}
  //}
}
