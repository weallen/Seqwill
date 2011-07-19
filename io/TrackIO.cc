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
  hid_t root_group = H5Gopen2(h5file_, "/", H5P_DEFAULT);
  hid_t trackgroup = H5Gopen2(root_group, trackname.c_str(), H5P_DEFAULT);
  H5G_info_t group_info;
  char* name;
  ssize_t size;

  if (H5Gget_info(trackgroup, &group_info) < 0) {
      WARNLOG("Can't find track " + trackname);
      return subtracknames;
  }
 for (hsize_t i = 0; i < group_info.nlinks; ++i) {
    // Have to call once with null name to get size
    size = H5Lget_name_by_idx(trackgroup, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, NULL, NULL, H5P_DEFAULT);
    name = new char[size+1];
    H5Lget_name_by_idx(trackgroup, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name, size+1, H5P_DEFAULT);
    std::string temp(name);
    subtracknames.push_back(temp);
    delete name;
  }
  H5Gclose(trackgroup);
  H5Gclose(root_group);
 return subtracknames;
}

std::vector<std::string>
TrackFile::GetTrackNames() const
{
  std::vector<std::string> tracknames;
  hid_t root_group = H5Gopen2(h5file_, "/", H5P_DEFAULT);
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

TrackMetadata
TrackFile::GetSubTrackMetadata(const std::string& trackname, const std::string& subtrackname) const
{
  TrackMetadata m;
  ScopedH5GOpen root_group(h5file_, "/");
  ScopedH5GOpen track_group(root_group.id(), trackname);
  if (track_group.id() < 0) {
    ERRORLOG("Couldn't load metadata for " + subtrackname);
    return m;
  }
  hid_t subtrack = H5Dopen2(track_group.id(), subtrackname.c_str(), H5P_DEFAULT);
  if (subtrack < 0) {
    ERRORLOG("Couldn't load metadata for subtrack " + subtrackname);
    return m;
  }

  ScopedH5AOpen start_attr(subtrack, "Start");
  ScopedH5AOpen stop_attr(subtrack, "Stop");
  ScopedH5AOpen res_attr(subtrack, "Resolution");
  int start;
  int stop;
  int resolution;
  H5Aread(start_attr.id(), H5T_NATIVE_INT, &start);
  H5Aread(stop_attr.id(), H5T_NATIVE_INT, &stop);
  H5Aread(res_attr.id(), H5T_NATIVE_INT, &resolution);
  m.SetIntMetadata("Start", start);
  m.SetIntMetadata("Stop", stop);
  m.SetIntMetadata("Resolution", resolution);
  H5Dclose(subtrack);
  return m;
}
