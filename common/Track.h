#ifndef TRACK_H_
#define TRACK_H_

#include <string>
#include <algorithm>
#include <boost/intrusive_ptr.hpp>
#include <Eigen/Core>

#include "data/TrackData.h"
#include "base/RefBase.h"
#include "base/Common.h"
#include "common/TrackMetadata.h"

template<typename DataT> class Track;

// Set of tracks attached to one genomic interval
template<typename DataType>
class Track : public RefBase
{
public:

  typedef DataType DataT;
  typedef boost::intrusive_ptr<Track<DataT> > Ptr;
  typedef std::vector<DataT> VectorType;

  Track() : RefBase() {}

  virtual ~Track() {}

  void Init(const std::string& name, int start, int stop, 
            const TrackMetadata& metadata) 
  {
      name_ = name;
      start_ = start;
      stop_ = stop;
      metadata_ = metadata;
      SizeChanged();
  }

  void SetData(const VectorType& td)
  {
    data_.assign(td.begin(), td.end());
  }

  // Iterator
  inline DataT operator[](size_t pos) { return data_[pos]; }
  inline size_t size() { return (data_.size()); }
  inline VectorType::iterator begin() { return data_.begin(); }
  // TODO Add interval iterator
  inline VectorType::iterator end() { return data_.end(); }

  void Clear(const DataT& val) 
  { std::fill(data_.begin(), data_.end(), val); }

  long long int MemSize() const 
  { return (sizeof(*this) + data_.capacity() * sizeof(DataT)); }

  // Getters and setters
  int start() const { return start_; }
  int stop() const { return stop_; }
  const std::string& name() const { return name_; }
  const VectorType& data() const { return data_; }
  const TrackMetadata& metadata() const { return metadata_; }

  void set_extends(int start, int stop) 
  { stop_ = stop; start_ = start; SizeChanged(); }
  void set_name(const std::string& name) { name_ = name; }
  void set_data(const VectorType& data) { data_ = data; }
  void set_metadata(const TrackMetadata& metadata) { metadata_ = metadata; }

private:
  DISALLOW_COPY_AND_ASSIGN(Track)

  void SizeChanged()
  {
      if (start_ < stop_) {
          data_.resize(stop_ - start_);
      }
  }

  int start_;
  int stop_;
  std::string name_;
  std::string chr_;
  VectorType data_;
  TrackMetadata metadata_;
};

template <typename DataT>
std::ostream& operator << (std::ostream& s, const Track<DataT>& t) {
  s << "Track: " << t.GetName() << " chr " <<
       t.chr() << " : " << t.GetStart() << " - " << t.GetEnd();
   
}

#endif
