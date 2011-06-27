#ifndef TRACKDATA_H_
#define TRACKDATA_H_

#include <string>
#include <vector>
#include <map>

// Idea: Template TrackData some data can be more complex that just float
// E.g. store +/- strand information

template <typename DataT>
class TrackData
{
public:

  std::string track_name;
  std::vector<DataT> track;

  std::string chr_name;
  int start;
  int end;

  TrackData()
  {}

  virtual ~TrackData() {}
};
#endif 
