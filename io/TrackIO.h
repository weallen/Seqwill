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

#include <hdf5.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "base/Common.h"
#include "base/RefCount.h"
#include "base/Hdf5Util.h"
#include "common/Track.h"
#include "io/Traits.h"

// A trackfile can store multiple subtracks
// For example, in genomic data these subtracks would correspond to individual chromosomes

class TrackIO : public RefBase
{
public:

  typedef boost::intrusive_ptr<TrackIO> Ptr;

  TrackIO() 
    : RefBase()
    , h5file_(-1) 
  {}

  virtual ~TrackIO() {
    Close();
  }

  bool Create(const char* fname);
  bool Create(const std::string& fname);
  bool Open(const char* fname);
  bool Open(const std::string& fname);
  void Close();
  bool IsOpen() { return isopen_; }

  // General Track stuff
  std::vector<std::string> GetSubTrackNames(const std::string& subtrackname) const;
  std::vector<std::string> GetTrackNames() const;

  // Write Stuff
  template <typename T>
  bool WriteSubTrack(const std::string& tname, typename Track<T>::Ptr data);
 
  // Read stuff
  template <typename T>
  bool ReadSubTrack(const std::string& trackname,
                    const std::string& subtrackname,
                    typename Track<T>::Ptr subtrack); 
private:
  
  hid_t h5file_;
  std::string filename_;
  bool isopen_;
};

// Implementation
// FUCK YOU C++



template<typename DataT>
bool TrackIO::WriteSubTrack(const std::string& trackname,
                            typename Track<DataT>::Ptr subtrack)
{
  hid_t dataset;
  hid_t track_group;
  hid_t root_group;
  std::vector<std::string> tracknames;
  H5SCreate data_space(H5S_SIMPLE);
  hid_t dcpl;

  if (data_space.id() < 0) {
    ERRORLOG("Can't create data space");
    return false;
  }


  hsize_t memsize[1] = {(hsize_t) subtrack->MemSize()};
  hsize_t best_chunk_size = 4096*16;
  const hsize_t chunk_size = std::min(best_chunk_size, memsize[0] / 2);

  root_group = H5Gopen2(h5file_, "/", H5P_DEFAULT);

  // Make sure track exists.
  tracknames = GetTrackNames();
  if (std::find(tracknames.begin(), tracknames.end(), trackname) == tracknames.end()) {
    track_group = H5Gcreate2(root_group, trackname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    track_group = H5Gopen2(root_group, trackname.c_str(), H5P_DEFAULT);
  }
  if (track_group < 0) {
    ERRORLOG("Can't open track group " + trackname);
    H5Gclose(root_group);
    return false;
  }

  H5Sset_extent_simple(data_space.id(), 1, memsize, NULL);
  // Check if subtrack exists already
  tracknames = GetSubTrackNames(trackname);
  if (std::find(tracknames.begin(), tracknames.end(), subtrack->name()) == tracknames.end()) {
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    dataset = H5Dcreate(track_group, subtrack->name().c_str(),
                        DataTypeTraits<DataT>::H5Type(), H5P_DEFAULT, dcpl);
  } else {
    dataset = H5Dopen2(track_group, subtrack->name().c_str(), H5P_DEFAULT);
  }

  if (dataset < 0) {
    ERRORLOG("Couldn't open subtrack " + subtrack->name());
    H5Gclose(root_group);
    H5Gclose(track_group);
    return false;
  }

  bool attrwritesuccess = true;
  // Write size attributes
  if (!WriteAttribute(dataset, "Start", 1, subtrack->start())) {
    attrwritesuccess = false;
  }

  if (!WriteAttribute(dataset, "Stop", 1, subtrack->stop())) {
    attrwritesuccess = false;
  }

  if (!WriteAttribute(dataset, "Name", subtrack->name())) {
    attrwritesuccess = false;
  }

  if (!attrwritesuccess) {
    ERRORLOG("Couldn't write attributes");
    H5Gclose(root_group);
    H5Gclose(track_group);
    H5Dclose(dataset);
    return false;
  }

  hid_t err = H5Dwrite(dataset, DataTypeTraits<DataT>::H5Type(), H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, &(*subtrack->begin()));
  if (err < 0) {
    ERRORLOG("Error writing subtrack");
    H5Gclose(root_group);
    H5Dclose(dataset);
    H5Gclose(track_group);
    return false;
  }
  H5Gclose(root_group);
  H5Dclose(dataset);
  H5Gclose(track_group);
  return true;
}
template <typename DataT>
bool
TrackIO::ReadSubTrack(const std::string& trackname,
                             const std::string& subtrackname,
                             typename Track<DataT>::Ptr subtrack)
{
  hid_t root_group;
  hid_t track_group;
  hid_t dataset;
  int start;
  int stop;

  root_group = H5Gopen2(h5file_, "/", H5P_DEFAULT);
  // Do some error checking...
  track_group = H5Gopen2(root_group, trackname.c_str(), H5P_DEFAULT);
  if (track_group < 0) {
    ERRORLOG("Can't open track group " + trackname);
    H5Gclose(root_group);
    return false;
  }
  dataset = H5Dopen2(track_group, subtrackname.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    ERRORLOG("Can't open subtrack " + subtrackname);
    H5Gclose(root_group);
    H5Gclose(track_group);
    return false;
  }

  // Read size attributes...
  if (!ReadAttribute(dataset, "Start", 1, &start)
      || !ReadAttribute(dataset, "Stop", 1, &stop)) {
    ERRORLOG("Can't read attributes");
    H5Gclose(root_group);
    H5Gclose(track_group);
    H5Dclose(dataset);
    return false;
  }
  subtrack->set_extends(start, stop);
  subtrack->set_name(subtrackname);
  if (H5Dread(dataset, DataTypeTraits<DataT>::H5Type(), H5S_ALL, H5S_ALL,
              H5P_DEFAULT, &(*subtrack->begin())) < 0) {
    ERRORLOG("Can't read dataset " + subtrackname);
    H5Gclose(root_group);
    H5Gclose(track_group);
    H5Dclose(dataset);
    return false;
  }
  H5Gclose(root_group);
  H5Gclose(track_group);
  H5Dclose(dataset);
  return true;
}


    
#endif
