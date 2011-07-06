#ifndef TRACK_METADATA_H_
#define TRACK_METADATA_H_

#include <string>
#include <map>

#include "base/Common.h"
#include "base/Types.h"

class TrackMetadata
{
public:
  typedef std::map<std::string, std::string> StrMetadata;
  typedef std::map<std::string, int> IntMetadata;
  typedef std::map<std::string, float> FloatMetadata;

  TrackMetadata() {}
  virtual ~TrackMetadata() {}

  float GetFloatMetadata(const std::string& name, const float defaultval) const
  {
    float ret = defaultval;
    FloatMetadata::const_iterator i = float_metadata_.find(name);
    if (i != float_metadata_.end()) {
        ret = i->second;
    }
    return ret;
  }

  int GetIntMetadata(const std::string& name, const int defaultval) const
  {
    int ret = defaultval;
    IntMetadata::const_iterator i = int_metadata_.find(name);
    if (i != int_metadata_.end()) {
        ret = i->second;
    }
    return ret;
  }


  std::string GetStrMetadata(const std::string& name, 
                             const std::string& defaultval) const
   {
    std::string ret = defaultval;
    StrMetadata::const_iterator i = str_metadata_.find(name);
    if (i != str_metadata_.end()) {
        ret = i->second;
    }
    return ret;
  }


  void SetFloatMetadata(const std::string& name, const float val)
  {
    float_metadata_[name] = val;
  }
   
  void SetIntMetadata(const std::string& name, const int val)
  {
    int_metadata_[name] = val;
  }

  void SetStrMetadata(const std::string& name, const std::string& val)
  {
    str_metadata_[name] = val;
  }

  
  const FloatMetadata& float_metadata() const 
  { return float_metadata_; }

  const StrMetadata& str_metadata() const
  { return str_metadata_; }

  const IntMetadata& int_metadata() const
  { return int_metadata_; }
 

private:
  DISALLOW_COPY_AND_ASSIGN(TrackMetadata)

  StrMetadata str_metadata_;
  IntMetadata int_metadata_;
  FloatMetadata float_metadata_;
};
#endif
