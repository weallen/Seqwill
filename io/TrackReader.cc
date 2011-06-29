#include "io/TrackReader.h"

void TrackReader::Open(const char* dirname)
{
  std::string s(dirname);
  Open(s);
}

void TrackReader::Open(const std::string& dirname)
{
}

void TrackReader::Close()
{
  Flush();
}

// Deallocate everything
void TrackReader::Flush()
{
}


// Chromosome stuff
std::vector<string> TrackReader::GetTrackNames()
{
  if (!isopen_) {
    Open();
  }
  return chrnames_;
}

bool TrackReader::HasChromosome(const string& chrname) {
  for (svec<string>::iterator i = chrnames_.begin();
      i != chrnames_.end(); ++i) {
    if (chrname == *i)
      return true;
  }
  return false;
}

TrackPtr TrackReader::GetChromosome(const std::string& chrname, const std::string& trackname)
{
  TrackPtr track(new Track(trackname, chrname, 0, chrlens_[chrname]));
  TrackData td;
  track_reader_.GetChromosome();
}

TrackPtr TrackReader::GetChrRegion(const std::string& chrname, int start, int end, const std::string& trackname)
{

}
