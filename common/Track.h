#ifndef TRACK_H_
#define TRACK_H_

#include <string>
#include <algorithm>
#include <math.h>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <Eigen/StdVector>

#include "data/TrackData.h"
#include "base/Common.h"
#include "common/TrackMetadata.h"
#include "base/Log.h"
#include "io/Traits.h"

template<typename DataT> class Track;

//
// Provide Read Only access to data.
template<typename T>
class Track
{
public:

  typedef T DataT;
  typedef boost::shared_ptr<Track<DataT> > Ptr;
  typedef boost::shared_ptr<const Track<DataT> > ConstPtr;
  typedef std::vector<DataT, Eigen::aligned_allocator<DataT> > VectorType;
  typedef typename VectorType::const_iterator const_iterator;
  typedef typename VectorType::iterator iterator;

  Track()
    : resolution_(1)
    , start_(-1)
    , stop_(-1)
    , track_("")
    , subtrack_("")
  {}

  Track(int start, int stop, const VectorType& data)
    : resolution_(1)
    , start_(start)
    , stop_(stop)
    , track_("")
    , subtrack_("")
    , data_(data)
  {}

  virtual ~Track() {}

  void Init(const std::string& name, int start, int stop, 
            const TrackMetadata& metadata) 
  {
    resolution_ = 1;
    track_ = name;
    start_ = start;
    stop_ = stop;
    data_.resize(stop_ - start_);
  }

  // Iterator

  inline DataT operator[](size_t pos) { return data_[pos]; }
  inline size_t size() { return (data_.size()); }
  inline const_iterator cbegin() const { return (data_.begin()); }
  inline const_iterator cend() const { return (data_.end()); }
  inline DataT& cvalue(size_t pos) const { return data_[pos]; }
  inline size_t size() const { return (data_.size()); }
  inline iterator begin() { return data_.begin(); }
  inline iterator end() { return data_.end(); }

  void Clear(const DataT& val)
  { std::fill(data_.begin(), data_.end(), val); }

  long long int MemSize() const 
  { return (sizeof(*this) + data_.capacity() * sizeof(DataT)); }

  // Getters and setters

  int start() const { return start_; }
  int stop() const { return stop_; }
  const std::string& trackname() const { return track_; }
  const std::string& subtrackname() const { return subtrack_; }
  const TrackMetadata& metadata() const { return metadata_; }
  int resolution() const { return resolution_; }

  void set_trackname(const std::string& name) { track_ = name; }
  void set_subtrackname(const std::string& name) { subtrack_ = name; }
  void set_metadata(const TrackMetadata& metadata) { metadata_ = metadata; }
  void set_resolution(int res) { resolution_ = res; }
  void set_data(const typename Track<DataT>::VectorType& data) { data_ = data; }
  void ReduceResolution(int rez);

  // these don't do anything unless it is a resizabletrack
  void set_extends(int start, int stop)
  { stop_ = stop; start_ = start; SizeChanged(); }

  virtual std::string ClassName() const
  { return std::string("Track"); }

  virtual std::string DataType() const
  { return DataTypeTraits<DataT>::Name(); }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
private:
  void SizeChanged()
  {
      if (start_ < stop_) {
          data_.resize(stop_ - start_);
      }
  }

  int ComputeResizedLen(int size)
  {
    int full_len = static_cast<float>(data_.size()) * resolution_;
    int new_len = ceil(full_len / size);
    return new_len;
  }

protected:
  DISALLOW_COPY_AND_ASSIGN(Track)


  int resolution_;
  int start_;
  int stop_;
  std::string track_;
  std::string subtrack_;
  typename Track<DataT>::VectorType data_;
  TrackMetadata metadata_;

};

template <typename DataT>
std::ostream& operator << (std::ostream& s, const Track<DataT>& t) {
  s << "Track: " << t.GetName() << " chr " <<
       t.chr() << " : " << t.GetStart() << " - " << t.GetEnd();
    return s;
}

template <typename T>
void Track<T>::ReduceResolution(int size)
{
  if ((size > resolution_) && (DataType().compare("char") != 0)) {
    resolution_ = size;
    int new_len = ComputeResizedLen(size);
    typename Track<T>::VectorType v(new_len);
    for (int i = 0; i < new_len; ++i) {
      v[i] = 0;
      for (int j = 0; j < size; ++j) {
        v[i] += data_[i * size + j];
      }
    }
    start_ = ceil(start_ / size);
    stop_ = ceil(stop_ / size);
    data_ = v;
  } else {
    WARNLOG("New resolution must be less than old resolution");
  }
}

#endif
