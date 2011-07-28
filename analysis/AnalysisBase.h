#ifndef ANALYSIS_BASE_H_
#define ANALYSIS_BASE_H_
#include <string>
#include <Eigen/StdVector>
#include <boost/shared_ptr.hpp>

#include "common/Track.h"

template <typename TypeT>
class AnalysisBase
{
public:
  typedef Track<TypeT> Track;
  typedef typename Track::Ptr TrackPtr;
  typedef typename Track::ConstPtr TrackConstPtr;

  AnalysisBase()
    : input_()
  {}

  virtual ~AnalysisBase();

  virtual void set_input_track(const TrackConstPtr& track) { input_ = track; }
  TrackConstPtr const  input_track() { return input_; }
  
  virtual void Compute() = 0;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
  bool Init()
  {
    if (!input_) {
      return false;
    }
    return true;
  }

  bool DeInit()
  {
    return true;
  }

  TrackPtr input_;
};

template <typename TrackInT, typename TrackOutT>
class Analysis : public AnalysisBase<TrackInT>
{
using AnalysisBase<TrackInT>::Init;
using AnalysisBase<TrackInT>::DeInit;

public:

  typedef AnalysisBase<TrackInT> BaseClass;
  typedef boost::shared_ptr<Analysis<TrackInT, TrackOutT> > Ptr;
  typedef boost::shared_ptr<const Analysis<TrackInT, TrackOutT> > ConstPTr;

  typedef Track<TrackInT> TrackIn;
  typedef typename Track<TrackInT>::Ptr TrackInPtr;
  typedef typename Track<TrackInT>::ConstPtr TrackInConstPtr;

  typedef Track<TrackOutT> TrackOut;
  typedef typename Track<TrackOutT>::Ptr TrackOutPtr;

  Analysis() 
  {}

  virtual ~Analysis() {}

  void Compute(TrackOutPtr output) 
  { ComputeAnalysis(output); }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
  // Set this is the leaf classes
  std::string analysis_name_;

  inline const std::string&
  ClassName() const { return analysis_name_; }

private:
  virtual void ComputeAnalysis(TrackOutPtr output) = 0;
};


template <typename TrackInT>
class Filter : public AnalysisBase<TrackInT>
{
public:
  using AnalysisBase<TrackInT>::Init;
  using AnalysisBase<TrackInT>::DeInit;

  typedef AnalysisBase<TrackInT> BaseClass;
  typedef boost::shared_ptr<Filter<TrackInT> > Ptr;
  typedef boost::shared_ptr<const Filter<TrackInT> > ConstPtr;

  typedef typename Track<TrackInT>::Ptr TrackPtr;
  typedef typename Track<TrackInT>::ConstPTr TrackConstPtr;

  Filter() {}
  virtual ~Filter() {}

  // Calls filter method 
  virtual void 
  Compute(typename Track<TrackInT>::Ptr out)
  {
    if (!Init())
      return;
    ApplyFilter(out);
    DeInit();
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
  virtual void
  ApplyFilter(typename Track<TrackInT>::Ptr out) = 0;

  inline const std::string& ClassName() const { return filter_name_; }
  std::string filter_name_;
};
#endif
