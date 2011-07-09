#include "io/TrackIO.h"
bool TrackFile::HasTrack(const std::string& trackname) const {
  std::vector<std::string> tracknames = GetTrackNames();
  if (std::find(tracknames.begin(), tracknames.end(), trackname) != tracknames.end()) {
    return true;
  }
  return false;
}

bool TrackFile::HasSubTrack(const std::string& trackname,
                          const std::string& subtrackname) const {
  std::vector<std::string> tracknames = GetSubTrackNames(trackname);
  if(std::find(tracknames.begin(), tracknames.end(), subtrackname) != tracknames.end()) {
    return true;
  }
  return false;
}

bool TrackFile::Create(const char* fname)
{
  return Create(std::string(fname));
}

bool TrackFile::Create(const std::string& fname)
{
  h5file_ = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h5file_ < 0) {
    ERRORLOG("Couldn't create file " + fname);
    return false;
  }
  isopen_ = true;
  return true;
}

bool TrackFile::Open(const char* dirname)
{
  return Open(std::string(dirname));
}

bool TrackFile::Open(const std::string& dirname)
{
  filename_ = dirname;
  bool success = true;
  struct stat s;
  if (stat(dirname.c_str(), &s) == 0) {
    if (H5Fis_hdf5(dirname.c_str())) {
      h5file_ = H5Fopen(dirname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      if (h5file_ < 0)
        return false;
    } else {
      DEBUGLOG("File " + dirname + " is not an HDF5 file");
      return false;
    }
    isopen_ = true;
  } else {
    ERRORLOG("File " + dirname + " does not exist.");
    return false;
  }
  return true;
}

void TrackFile::Close()
{
  if (isopen_) {
    isopen_ = false;
    H5Fclose(h5file_);
  }
}

std::vector<std::string>
TrackFile::GetSubTrackNames(const std::string& trackname) const
{
  std::vector<std::string> subtracknames;
  ScopedH5GOpen trackgroup(h5file_, trackname);
  hsize_t idx;
  H5G_info_t group_info;
  char* name;
  ssize_t size;

  if (H5Gget_info(trackgroup.id(), &group_info) < 0) {
      WARNLOG("Can't find track " + trackname);
      return subtracknames;
  }
 for (hsize_t i = 0; i < group_info.nlinks; ++i) {
    // Have to call once with null name to get size
    size = H5Lget_name_by_idx(trackgroup.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, NULL, NULL, H5P_DEFAULT);
    name = new char[size+1];
    H5Lget_name_by_idx(trackgroup.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name, size+1, H5P_DEFAULT);
    std::string temp(name);
    subtracknames.push_back(temp);
    delete name;
  }
 return subtracknames;
}

std::vector<std::string>
TrackFile::GetTrackNames() const
{
  std::vector<std::string> tracknames;
  hid_t root_group = H5Gopen2(h5file_, "/", H5P_DEFAULT);
  hsize_t idx;
  H5G_info_t group_info;
  char* name;
  ssize_t size;

  if (H5Gget_info(root_group, &group_info) < 0) {
      WARNLOG("Can't find root group");
      return tracknames;
  }
 for (hsize_t i = 0; i < group_info.nlinks; ++i) {
    // Have to call once with null name to get size
    size = H5Lget_name_by_idx(root_group, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, NULL, NULL, H5P_DEFAULT);
    name = new char[size+1];
    H5Lget_name_by_idx(root_group, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name, size+1, H5P_DEFAULT);
    std::string temp(name);
    tracknames.push_back(temp);
    delete name;
  }
 H5Gclose(root_group);
 return tracknames;
}
