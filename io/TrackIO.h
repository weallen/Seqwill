#ifndef TRACKIO_H_
#define TRACKIO_H_
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <boost/intrusive_ptr.hpp>
#include <boost/filesystem.hpp>

#include <hdf5.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "common/Track.h"
#include "base/Common.h"
#include "base/RefBase.h"

template <typename TypeT>
int LoadTrack(const std::string& fname, Track<TypeT>* track);


template <typename T>
int SaveTrack(const std::string& fname, Track<T>* track);

// A trackfile can store multiple subtracks
// For example, in genomic data these subtracks would correspond to individual chromosomes
template<typename T>
class TrackIO : public RefBase
{
public:

  typedef boost::intrusive_ptr<TrackIO<T> > Ptr;

  TrackIO() 
    : RefBase()
    , h5file_(-1) 
  {}

  virtual ~TrackIO() {
    Close();
  }

  bool Open(const char* fname);
  bool Open(const std::string& fname);
  void Close();
  bool IsOpen() { return isopen_; }
  void Flush();

  // General Track stuff
  std::vector<std::string> GetSubTrackNames(const std::string& subtrackname) const;
  std::vector<std::string> GetTrackNames() const;

  // Write Stuff
  bool WriteSubTrack(const std::string& tname, typename Track<T>::Ptr data);
 
  // Read stuff
  
  typename Track<T>::Ptr ReadSubTrack(const std::string& trackname,
                                      const std::string& subtrackname); 
private:
  void LoadTrackNames();

  hid_t h5file_;
  std::string filename_;
  bool isopen_;
};

    
#endif
