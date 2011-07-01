#ifndef TRACKREADER_H_
#define TRACKREADER_H_
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include <hdf5.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "common/Track.h"
#include "hdf5/HDFTrackReader.h"
#include "hdf5/HDFTrackWriter.h"

template <typename TypeT>
int LoadTrack(const std::string& fname, Track<TypeT>* track);

// A trackfile can store multiple subtracks
// For example, in genomic data these subtracks would correspond to individual chromosomes
template<typename T>
class TrackReader
{
public:

  TrackReader() {}

  virtual ~TrackReader() {
    Close();
  }

  void Open(const char* fname);
  void Open(const std::string& fname);
  void Close();
  bool IsOpen() { return isopen_; }
  void Flush();


  // Track stuff
  std::vector<std::string> GetSubTrackNames();
  std::string GetTrackName();

  // Caller is responsible for allocating and deallocating memory for tracks.
  int ReadSubTrack(const std::string& subtrackname, Track<T>* track);

  // Return -1 on error, 1 on success
  int ReadSubTrackRegion(const std::string& subtrackname, int start, int end, Track<T>* track);



private:
  void LoadTrackNames();

  string filename_;
  bool isopen_;

  HDFTrackReader<T> track_reader_;
  std::vector<string> chrnames_;
  std::map<std::string, int> chrlens_;
  std::vector<string> tracknames_;
};

#endif
