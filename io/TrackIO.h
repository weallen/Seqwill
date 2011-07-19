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

#include <boost/shared_ptr.hpp>

#include <hdf5.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "base/Common.h"
#include "base/RefCount.h"
#include "base/Hdf5Util.h"
#include "common/Track.h"
#include "io/Traits.h"
#include "common/TrackMetadata.h"


class TrackFile;


// A trackfile can store multiple subtracks
// For example, in genomic data these subtracks would correspond to individual chromosomes
class TrackFile
{
public:

  typedef boost::shared_ptr<TrackFile> Ptr;

  TrackFile()
    : h5file_(-1)
  {}
   
  virtual ~TrackFile() {
    Close();
  }


  // General Track stuff
  std::vector<std::string> GetSubTrackNames(const std::string& subtrackname) const;
  std::vector<std::string> GetTrackNames() const;
  TrackMetadata GetSubTrackMetadata(const std::string& trackname, const std::string& subtrackname) const;
  
  bool HasSubTrack(const std::string& trackname, const std::string& subtrackname) const;
  bool HasTrack(const std::string& trackname) const;

  // Write Stuff
  template <typename T>
  bool WriteSubTrack(const std::string& fname,
                     const std::string& tname,
                     typename Track<T>::Ptr data);
 
  // Read stuff
  template <typename T>
  bool ReadSubTrack(const std::string& fname,
                    const std::string& trackname,
                    const std::string& subtrackname,
                    typename Track<T>::Ptr subtrack);
  
  bool Open(const char* fname);
  bool Open(const std::string& fname);
  void Close();
  bool IsOpen() { return isopen_; }

private:
  DISALLOW_COPY_AND_ASSIGN(TrackFile)

  bool Create(const char* fname);
  bool Create(const std::string& fname);

  hid_t h5file_;
  std::string filename_;
  bool isopen_;
};

// Implementation
// FUCK YOU C++

//-------------------------------------------------------------------
// Implementation

template<typename DataT>
bool TrackFile::WriteSubTrack(const std::string& fname,
                            const std::string& trackname,
                            typename Track<DataT>::Ptr subtrack)
{
  hid_t dataset;
  hid_t track_group;
  hid_t root_group;
  std::vector<std::string> tracknames;
  ScopedH5SCreate data_space(H5S_SIMPLE);
  hid_t dcpl;

  if (!Open(fname)) {
    if (!Create(fname)) {
      ERRORLOG("Can't open file " + fname);
      return false;
    }
  }

  if (data_space.id() < 0) {
    ERRORLOG("Can't create data space");
    return false;
  }


  hsize_t memsize[1] = {(hsize_t) subtrack->MemSize()};
  //hsize_t best_chunk_size = 4096*16;
//  const hsize_t chunk_size = std::min(best_chunk_size, memsize[0] / 2);

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
  if (std::find(tracknames.begin(), tracknames.end(), subtrack->subtrackname()) == tracknames.end()) {
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    dataset = H5Dcreate2(track_group, subtrack->subtrackname().c_str(),
                        DataTypeTraits<DataT>::H5Type(), data_space.id(),
                        H5P_DEFAULT, dcpl, H5P_DEFAULT);
  } else {
    dataset = H5Dopen2(track_group, subtrack->subtrackname().c_str(), H5P_DEFAULT);
  }

  if (dataset < 0) {
    ERRORLOG("Couldn't open subtrack " + subtrack->subtrackname());
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

  if (!WriteAttribute(dataset, "Name", subtrack->subtrackname())) {
    attrwritesuccess = false;
  }

  if (!WriteAttribute(dataset, "Resolution", 1, subtrack->resolution())) {
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
  Close();
  return true;
}
template <typename DataT>
bool
TrackFile::ReadSubTrack(const std::string& fname,
                      const std::string& trackname,
                      const std::string& subtrackname,
                      typename Track<DataT>::Ptr subtrack)
{
  hid_t dataset;
  int start;
  int stop;
  int resolution;

  if(!Open(fname)) {
    ERRORLOG("Can't open " + fname);
    return false;
  }
  ScopedH5GOpen root_group(h5file_, std::string("/"));
  // Do some error checking...
  ScopedH5GOpen track_group(root_group.id(), trackname);
      //= H5Gopen2(root_group, trackname.c_str(), H5P_DEFAULT);
  if (track_group < 0) {
    ERRORLOG("Can't open track group " + trackname);
    return false;
  }
  // CHeck if dataset exists before trying to open
  if (!HasSubTrack(trackname, subtrackname)) {
    WARNLOG("Can't find subtrack " + subtrackname);
    return false;
  }
  dataset = H5Dopen2(track_group, subtrackname.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    ERRORLOG("Can't open subtrack " + subtrackname);
    return false;
  }

  // Read size attributes...
  if (!ReadAttribute(dataset, "Start", 1, &start)
      || !ReadAttribute(dataset, "Stop", 1, &stop)
      || !ReadAttribute(dataset, "Resolution", 1, &resolution)) {
    ERRORLOG("Can't read attributes");
    H5Dclose(dataset);
    return false;
  }
  subtrack->set_resolution(resolution);
  subtrack->set_extends(start, stop);
  subtrack->set_trackname(trackname);
  subtrack->set_subtrackname(subtrackname);
  if (H5Dread(dataset, DataTypeTraits<DataT>::H5Type(), H5S_ALL, H5S_ALL,
              H5P_DEFAULT, &(*subtrack->begin())) < 0) {
    ERRORLOG("Can't read dataset " + subtrackname);
    H5Dclose(dataset);
    return false;
  }
  H5Dclose(dataset);
  return true;
}


    
#endif
