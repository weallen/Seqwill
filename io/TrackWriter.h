#ifndef TRACKWRITER_H_
#define TRACKWRITER_H_

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
#include "common/Track.h"
#include "hdf5/HDFTrackReader.h"
#include "hdf5/HDFTrackWriter.h"

template <typename TypeT>
int SaveTrack(const std::string& fname, Track<TypeT>* track);

template <typename DataT>
class TrackWriter
{
public:
  TrackWriter();
  virtual ~TrackWriter()
  {
    Close();
  }

  void Open(const char* fname);
  void Open(const std::string& fname);



 // Sequence stuff
  void LoadChrSeq(const std::string& seqname);

  void LoadTrackData(const std::string& tdfile);
};
#endif
