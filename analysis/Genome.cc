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
GenomeData::GetSubTrackForChrom(const std::string& chrname)
{
  Track<float>::Ptr curr_track;

  if (open_chrs_.find(chrname) != open_chrs_.end()) {
    Track<float>::Ptr curr_track = open_chrs_[chrname];    
  } else {
    Track<float>::Ptr curr_track(new Track<float>);    
    trackfile_->ReadSubTrack<float>(trackname_, chrname, *curr_track);
    open_chrs_[chrname] = curr_track;
  }
  return curr_track;
}

void
GenomeData::SaveTrackFromWIG(const std::string& wigname, int resolution)
{
  WIGParser p;
  DEBUGLOG("Loading track " + trackname_ + " from WIG file " + wigname);    
  p.Open(wigname);
  std::vector<std::string> chrnames = genome_info_.chr_names();

  WIGLine curr_line;
  Track<float>::Ptr curr_track;
  std::string curr_chr;
  int max_pos;
  while(p.HasNextLine()) {
    curr_line = p.NextLine();
    if (p.StateChanged()) {
      if (curr_track.get() != NULL) {
        DEBUGLOG("Writing " + curr_track->subtrackname());
        trackfile_->WriteSubTrack<float>(*curr_track);
      }
      curr_chr = ChrToString(p.curr_state().chr);
      DEBUGLOG("Processing " + curr_chr);
      curr_track = open_chrs_[curr_chr];
      curr_track = Track<float>::Ptr(new Track<float>);
      curr_track->set_resolution(resolution);
      curr_track->set_abs_extends(0, genome_info_.chr_size(curr_chr));
      curr_track->set_trackname(trackname_);
      curr_track->set_subtrackname(curr_chr);
      curr_track->set_extends(0, ceil(static_cast<float>(genome_info_.chr_size(curr_chr))/resolution));      
      max_pos = curr_track->astop();
    }
    if (curr_line.pos < max_pos) {
      curr_track->aset(curr_line.pos, curr_line.val);    
    }
  }
  // Write last track
  DEBUGLOG("Writing " + curr_track->subtrackname());
  trackfile_->WriteSubTrack<float>(*curr_track);
}
