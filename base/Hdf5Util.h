#ifndef HDF5UTIL_H_
#define HDF5UTIL_H_
#include <string>
#include <vector>
#include <hdf5.h>

#include "base/Log.h"

class H5Base
{
public:
  H5Base()
    : id_(-1)
  {}

  hid_t id() const 
  { return id_; }

  operator hid_t ()
  { return id_; }

protected:
  hid_t id_;
};

// GROUP STUFF
class H5GOpen : public H5Base
{
public:
  H5GOpen()
    : H5Base() 
  {}

  H5GOpen(hid_t loc, const std::string& name) 
  {
    id_ = H5Gopen2(loc, name.c_str(), H5P_DEFAULT);
    if (id_ < 0) {
        ERRORLOG("Couldn't open attribute " + name);
    }
  }

  ~H5GOpen() {
    if (id_ >= 0) {
      H5Gclose(id_);
    }
  }
};

class H5GCreate : public H5Base
{
public:
  H5GCreate(hid_t loc, const std::string& name) 
  {
    id_ = H5Gcreate2(loc, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (id_ < 0) {
        ERRORLOG("Couldn't create group " + name);
    }
  }

  ~H5GCreate() {
    if (id_ >= 0) {
      H5Gclose(id_);
    }
  }
};

// ATTRIBUTE STUFF
class H5AOpen : public H5Base
{
public:
  H5AOpen(hid_t loc, const std::string& name) 
  {
    id_ = H5Aopen_name(loc, name.c_str());
    if (id_ < 0) {
        ERRORLOG("Couldn't open attribute" + name);
    }
  }

  ~H5AOpen() {
    if (id_ >= 0) {
      H5Aclose(id_);
    }
  }
};


class H5AOpenIdx : public H5Base
{
public:
  H5AOpenIdx(hid_t loc, unsigned int idx)
  {
    id_ = H5Aopen_idx(loc, idx);
    if (id_ < 0) {
      ERRORLOG("Could open attributed with index " + idx);
    }
  }

  ~H5AOpenIdx()
  {
    if (id_ >= 0) {
      H5Aclose(id_);
    }
  }
};

class H5AGetSpace : public H5Base
{
public:
  H5AGetSpace(hid_t attrid) 
  {
    id_ = H5Aget_space(attrid);
    if (id_ < 0) {
      ERRORLOG("Couldn't get attr space");
    }
  }

  ~H5AGetSpace()
  {
    if (id_ >= 0) {
      H5Sclose(id_);
    } 
  }
};

class H5AGetType : public H5Base
{
public:
  H5AGetType(hid_t attrid) 
  {
    id_ = H5Aget_space(attrid);
    if (id_ < 0) {
      ERRORLOG("Couldn't get attr space");
    }
  }

  ~H5AGetType()
  {
    if (id_ >= 0) {
      H5Tclose(id_);
    } 
  }
};


// DATASET STUFF


class H5DOpen : public H5Base
{
public:
  H5DOpen() : H5Base() {}
  H5DOpen(hid_t parent_loc, const std::string& name, hid_t dopl_id)
  {
    id_ = H5Dopen2(parent_loc, name.c_str(), dopl_id);
  }
  ~H5DOpen() 
  {
    if (id_ >= 0) {
      H5Dclose(id_);
    }
  }
};

class H5DCreate : public H5Base
{
public:
  H5DCreate(hid_t parent, const std::string& name,
            hid_t dtype_id, hid_t space_id, hid_t lcpl_id,
            hid_t dcpl_id, hid_t dapl_id)
  {
    id_ = H5Dcreate2(parent, name.c_str(), dtype_id, space_id, lcpl_id,
                    dcpl_id, dapl_id);
  }

  ~H5DCreate()
  {
    if (id_ >= 0) {
      H5Dclose(id_);
    }
  }
};

class H5DGetSpace : public H5Base
{
public:
    H5DGetSpace() : H5Base() {}
    H5DGetSpace(hid_t dsetid) 
    { 
      id_ = H5Dget_space(dsetid);
    }
    ~H5DGetSpace()
    {
      if (id_ >= 0) {
        H5Sclose(id_);
      }
    }
};

class H5DGetType : public H5Base
{
public:
    H5DGetType() : H5Base() {}
    H5DGetType(hid_t dsetid)
    { 
      id_ = H5Dget_type(dsetid);
    }
    ~H5DGetType()
    {
      if (id_ >= 0) {
        H5Tclose(id_);
      }
    }
};


// DATASPACE STUFF
class H5SCreate : public H5Base
{
public:
  H5SCreate() : H5Base() {}
  H5SCreate(H5S_class_t c)
  {
    id_ = H5Screate(c);
  }

  ~H5SCreate()
  {
    if (id_ >= 0) {
      H5Sclose(id_);
    }
  }
};

bool ReadAttribute(hid_t loc, const std::string& attr_name, std::string* value);
bool ReadAttribute(hid_t loc, const std::string& attr_name,
                   unsigned int attr_size, int* value);
bool ReadAttribute(hid_t loc, const std::string& attr_name,
                   unsigned int attr_size, float* value);

bool WriteAttribute(hid_t loc, const std::string& attr_name,
                    const std::string& value);
bool WriteAttribute(hid_t loc, const std::string& attr_name,
                    unsigned int attr_size, int value);
bool WriteAttribute(hid_t loc, const std::string& attr_name,
                    unsigned int attr_size, float value);
#endif
