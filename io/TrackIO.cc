#include "io/TrackIO.h"

bool TrackIO::Open(const char* dirname)
{
  return Open(std::string(dirname));
}

bool TrackIO::Open(const std::string& dirname)
{
  filename_ = dirname;
  bool success = true;
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
    ERRORLOG("Couldn't open file " + dirname);
    return false;
  }
  isopen_ = true;
  return true;
}

void TrackIO::Close()
{
  if (isopen_) {
    isopen_ = false;
    H5Fclose(h5file_);
  }
}

std::vector<std::string>
TrackIO::GetSubTrackNames(const std::string& trackname) const
{
  std::vector<std::string> subtracknames;
  H5GOpen trackgroup(h5file_, trackname);
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
    H5Lget_name_by_idx(trackgroup.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name, size, H5P_DEFAULT);
    std::string temp(name);
    subtracknames.push_back(temp);
  }
 return subtracknames;
}

std::vector<std::string>
TrackIO::GetTrackNames() const
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
    H5Lget_name_by_idx(root_group, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name, size, H5P_DEFAULT);
    std::string temp(name);
    tracknames.push_back(temp);
  }
 H5Gclose(root_group);
 return tracknames;
}
