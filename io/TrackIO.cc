#include "io/TrackIO.h"

template<typename T>
void TrackIO<T>::Open(const char* dirname)
{
  std::string s(dirname);
  Open(s);
}

template<typename T>
bool TrackIO<T>::Open(const std::string& dirname)
{
  filename_ = dirname;
  success = true;
  if (boost::filesystem::exists(dirname)) {
    if (H5Fis_hdf5(dirname.c_str())) {
      h5file_ = H5Fopen(dirname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      if (h5file_ < 0) 
        success = false;
    } else { 
        success = false; 
    }
  } else {
      success = false;
  }
  if (!success) {
    std::cerr << "Couldn't open file " + dirname;
    return false;
  }
  isopen_ = true;
  return true;
}

template<typename T>
void TrackIO<T>::Close()
{
  if (isopen_) {
    isopen_ = false;
    H5Fclose(h5file_);
  }
}

std::vector<std::string> GetSubTrackNames(const std::string& trackname) const
{
  std::vector<std::string> subtracknames;
  H5GOpen trackgroup(h5file_, trackname);
  hsize_t idx;
  H5G_info_t group_info;
  char* name;
  ssize_t size;
 
  if (H5Gget_info(trackgroup(), &group_info) < 0) {
      WARNLOG << "Can't find track " << trackname;
      return subtracknames;
  }
 for (hsize_t i = 0; i < group_info.nlinks; ++i) {
    // Have to call once with null name to get size
    size = H5Lget_name_by_idx(trackgroup(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, NULL, NULL, H5P_DEFAULT);
    H5Lget_name_by_idx(trackgroup(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, name, size, H5P_DEFAULT);
    std::string temp(name);
    subtracknames.push_back(temp);
  }
 return subtracknames;
}

std::vector<std::string> GetTrackNames() const
{
  std::vector<std::string> tracknames;
  hid_t rootgroup = H5Gopen(h5file_, "/");
  hsize_t idx;
  H5G_info_t group_info;
  char* name;
  ssize_t size;
 
  if (H5Gget_info(rootgroup, &group_info) < 0) {
      WARNLOG << "Can't find root group"; 
      return tracknames;
  }
 for (hsize_t i = 0; i < group_info.nlinks; ++i) {
    // Have to call once with null name to get size
    size = H5Lget_name_by_idx(trackgroup(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, NULL, NULL, H5P_DEFAULT);
    H5Lget_name_by_idx(trackgroup(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, name, size, H5P_DEFAULT);
    std::string temp(name);
    tracknames.push_back(temp);
  }
 H5Gclose(rootgroup);
 return tracknames;
}


template<typename DataT>
bool TrackIO::WriteSubTrack(const std::string& tname, 
                            typename Track<DataT>::Ptr subtrack)
{
  hsize_t memsize[1] = {(hsize_t) subtrack->MemSize()};
  hsize_t best_chunk_size = 4096*16;
  const hsize_t chunk_size = std::min(best_chunk_size, memsize[0] / 2);

  // Check if subtrack exists already

  // If subtrack does exist, just write to the region that the subtrack covers

  // If subtrack doesn't exist...
  hid_t err = H5Dwrite(dataset, DataTypeTraits<DataT>::h5type(), H5S_ALL, H5S_ALL, 
                       H5P_DEfAULT, &(*subtrack->begin()));
  if (err < 0) {
      std::cerr << "Error writing subtrack in TrackIO::WriteSubTrack";
      return false;
  }
  return true;
}
template <typename DataT>
typename Track<DataT>::Ptr subtrack
ReadSubTrack(const std::string& trackname, const std::string& subtrackname)
{
  // Read size attributes...
   
  return ReadSubTrackRegion(subtrackname, 0, end);
}

