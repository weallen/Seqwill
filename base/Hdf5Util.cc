#include "base/Hdf5Util.h"

bool ReadAttribute(hid_t loc, const std::string& attr_name, std::string& value)
{
  hsize_t str_len;
  H5A_info_t attr_info;
  H5T_class_t type_class;

  if (H5Aexists(loc, attr_name.c_str()) < 1) {
    ERRORLOG("Couldn't find attribute " + attr_name);
    return false;
  }
  ScopedH5AOpen attr(loc, attr_name);
  ScopedH5AGetSpace attr_space(attr);
  ScopedH5AGetType attr_type(attr);

  if (H5Aget_info(attr, &attr_info) != 1) {
    ERRORLOG("Bad attribute info for attribute " + attr_name);
    return false;
  } else {
      str_len = attr_info.data_size;
  }
  type_class = H5Tget_class(attr_type);
  if (type_class != H5T_STRING) {
      ERRORLOG("Wrong attribute type for attribute " + attr_name);
      return false;
  }
  std::vector<char> temp_str(str_len + 1);
  if (H5Aread(attr, H5T_C_S1, &temp_str[0]) < 0) {
    ERRORLOG("Couldn't read attribute " + attr_name);
    return false;
  }
  value = std::string(&temp_str[0]);
  return true;
}

bool ReadAttribute(hid_t loc, const std::string& attr_name, 
                   unsigned int attr_size, int* value)
{
  H5T_class_t type_class;

  if (H5Aexists(loc, attr_name.c_str()) < 1) {
    ERRORLOG("Couldn't find attribute " + attr_name);
    return false;
  }
  ScopedH5AOpen attr(loc, attr_name);
  ScopedH5AGetSpace attr_space(attr);
  ScopedH5AGetType attr_type(attr);
  if (H5Sget_simple_extent_ndims(attr_space.id()) != 1) {
    ERRORLOG("Bad attribute rank for attribute " + attr_name);
    return false;
  }
  hsize_t dims[1];
  H5Sget_simple_extent_dims(attr_space.id(), dims, NULL);
  if (dims[0] != attr_size) {
      ERRORLOG("Wrong attribute size for attribute " + attr_name);
      return false;
  }
  type_class = H5Tget_class(attr_type.id());
  if (type_class != H5T_INTEGER) {
      ERRORLOG("Wrong attribute type for attribute " + attr_name);
      return false;
  }
  if (H5Aread(attr, H5T_NATIVE_INT, value) < 0) {
    ERRORLOG("Couldn't read attribute " + attr_name);
    return false;
  }
  return true;
}

bool ReadAttribute(hid_t loc, const std::string& attr_name, 
                   unsigned int attr_size, float* value)
{
  H5T_class_t type_class;

  if (H5Aexists(loc, attr_name.c_str()) < 1) {
    ERRORLOG("Couldn't find attribute " + attr_name);
    return false;
  }
  ScopedH5AOpen attr(loc, attr_name);
  ScopedH5AGetSpace attr_space(attr);
  ScopedH5AGetType attr_type(attr);
  if (H5Sget_simple_extent_ndims(attr_space) != 1) {
    ERRORLOG("Bad attribute rank for attribute " + attr_name);
    return false;
  }
  hsize_t dims[1];
  H5Sget_simple_extent_dims(attr_space, dims, NULL);
  if (dims[0] != attr_size) {
      ERRORLOG("Wrong attribute size for attribute " + attr_name);
      return false;
  }
  type_class = H5Tget_class(attr_type);
  if (type_class != H5T_FLOAT) {
      ERRORLOG("Wrong attribute type for attribute " + attr_name);
      return false;
  }
  if (H5Aread(attr, H5T_NATIVE_FLOAT, value) < 0) {
    ERRORLOG("Couldn't read attribute " + attr_name);
    return false;
  }
  return true;
}

bool WriteAttribute(hid_t loc, const std::string& attr_name,
                    const std::string& value)
{
  hid_t attr = -1;
  hid_t attr_space;
  hid_t attr_type;
  bool success = true;
  attr_space = H5Screate(H5S_SCALAR);
  if (attr_space == -1)
   success = false;
  attr_type = H5Tcopy(H5T_C_S1);
  if (attr_type == -1)
    success = false;
  if (value.size()) {
    if (success && (H5Tset_size(attr_type, value.size()) == -1))
      success = false;
  }
  if (success) {
    H5Tset_strpad(attr_type, H5T_STR_NULLTERM);
    if (H5Aexists_by_name(loc, ".", attr_name.c_str(), H5P_DEFAULT) == 0) {
      attr = H5Acreate2(loc, attr_name.c_str(), attr_type,
                        attr_space, H5P_DEFAULT, H5P_DEFAULT);

    } else {
      attr = H5Aopen_by_name(loc, ".", attr_name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    }
  }

  if (attr == -1) {
    ERRORLOG("Couldn't create attribute " + attr_name);
    success = false;
  }
  if (success && (H5Awrite(attr, attr_type, value.c_str()) == -1)) {
    ERRORLOG("error writing attribute: " + attr_name);
    success = false;
  }

  
  H5Aclose(attr);
  H5Tclose(attr_type);
  H5Sclose(attr_space);
  return success;
}

bool WriteAttribute(hid_t loc, const std::string& attr_name, 
                    unsigned int attr_size, int value)
{
  hid_t attr;
  hid_t attr_space;
  hsize_t dims[1];
  dims[0] = attr_size;

  attr_space = H5Screate(H5S_SIMPLE);
  if (attr_space == -1)
   return false;
  if (H5Sset_extent_simple(attr_space, 1, dims, NULL) < 0)
    return false;
  if (H5Aexists_by_name(loc, ".", attr_name.c_str(), H5P_DEFAULT) == 0) {
    attr = H5Acreate2(loc, attr_name.c_str(), H5T_NATIVE_INT,
                      attr_space, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    attr = H5Aopen_by_name(loc, ".", attr_name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
  }
  if (attr < 0) {
    ERRORLOG("Couldn't create attribute " + attr_name);
    H5Aclose(attr);
    H5Sclose(attr_space);
    return false;
  }
  if (H5Awrite(attr, H5T_NATIVE_INT, &value) < 0) {
    ERRORLOG("error writing attribute: " + attr_name);
    H5Aclose(attr);
    H5Sclose(attr_space);
    return false;
  }

  
  H5Aclose(attr);
  H5Sclose(attr_space);
  return true;
}

bool WriteAttribute(hid_t loc, const std::string& attr_name,
                    unsigned int attr_size, float value)
{
  hid_t attr;
  hid_t attr_space;
  hsize_t dims[1];
  dims[0] = attr_size;

  attr_space = H5Screate(H5S_SIMPLE);
  if (attr_space == -1)
   return false;
  if (H5Sset_extent_simple(attr_space, 1, dims, NULL) < 0)
    return false;

  if (H5Aexists_by_name(loc, ".", attr_name.c_str(), H5P_DEFAULT) == 0) {
    attr = H5Acreate2(loc, attr_name.c_str(), H5T_NATIVE_FLOAT,
                      attr_space, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    attr = H5Aopen_by_name(loc, ".", attr_name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
  }


  if (attr < 0) {
    ERRORLOG( "Couldn't create attribute " + attr_name);
    H5Aclose(attr);
    H5Sclose(attr_space);
    return false;
  }
  if (H5Awrite(attr, H5T_NATIVE_INT, &value) < 0) {
    ERRORLOG("error writing attribute: " + attr_name);
    H5Aclose(attr);
    H5Sclose(attr_space);
    return false;
  }

  
  H5Aclose(attr);
  H5Sclose(attr_space);
  return true;
}


