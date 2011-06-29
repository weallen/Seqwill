#ifndef TRACK_H_
#define TRACK_H_

#include <string>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <Eigen/Array>

#include "data/TrackData.h"


// Set of tracks attached to one genomic interval
template<typename DataT>
class Track
{
public:

  typedef boost::shared_ptr<Track<DataT> > Ptr;
  typedef std::vector<DataT, Eigen::aligned_allocator<float> > VectorType;
  Track() {}
  Track(const TrackData<DataT>& td)
  {
    start = td->start;
    end = td->end;
    name = td->track_name;
    data.assign(td->track.begin(), td->track.end());
    chr = td->chr_name;
  }

  Track(const std::string& tname, const std::string& tchr, int ts, int te)
  {
    start = ts;
    end = te;
    name = tname;
    chr = tchr;
    data.resize(start - end);
  }

  virtual ~Track() {}

  void SetData(const std::vector<DataT>& td)
  {
    data.assign(td.begin(), td.end());
  }

  // Iterator
  inline operator[](size_t pos) { return data[pos]; }
  inline size_t size() { return (data.size()); }
  inline VectorType::iterator begin() { return data.begin(); }
  inline VectorType::iterator end() { return data.end(); }

  int start;
  int end;
  std::string name;
  std::string chr;
  VectorType data;
};

template <typename DataT>
std::ostream& operator << (std::ostream& s, const Track<DataT>& t) {
  s << "Track: " << t.GetName() << " chr " <<
       t.chr() << " : " << t.GetStart() << " - " << t.GetEnd();
}

#endif
