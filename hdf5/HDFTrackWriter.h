#ifndef HDFTRACKWRITER_H_
#define HDFTRACKWRITER_H_

#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include <boost/scoped_ptr.hpp>

#include "hdf5/HDFGroup.h"
#include "hdf5/HDFFile.h"
#include "hdf5/HDFAtom.h"
#include "hdf5/BufferedHDF2DArray.h"
#include "hdf5/BufferedHDFArray.h"

#include "base/StringUtil.h"
#include "base/DNASequence.h"

#include "data/TrackData.h"

typedef boost::scoped_ptr<BufferedHDF2DArray> Buffered2DArrayPtr;
typedef boost::scoped_ptr<BufferedHDFArray> BufferedArrayPtr;

class HDFTrackWriter
{
  public:

  void Create(std::string& outfilename) {
    outfilename_ = outfilename;
    outfile_.Create(outfilename_);
    root_group_.Initialize(*outfile_.hdfFile, "/");
    root_group_.AddGroup("chr");
    root_group_.AddGroup("seq");

    chr_group_.Initialize(root_group_.group, "chr");
    seq_group_.Initialize(root_group_.group, "seq");

    std::vector<std::string> empty_string;
    chr_lens_.Create(root_group_, "ChrLens");
    chr_lens_.Initialize(root_group_, "ChrLens");
    chr_lens_.Write(empty_string);
    chr_names_.Create(root_group_, "ChrNames");
    chr_names_.Initialize(root_group_, "ChrNames");
    chr_names_.Write(empty_string);
    track_names_.Create(root_group_, "TrackNames");
    track_names_.Initialize(root_group_, "TrackNames");
    track_names_.Write(empty_string);
  }

  void Open(std::string& outfilename)
  {
    outfilename_ = outfilename;
    outfile_.Create(outfilename_);
    root_group_.Initialize(*outfile_.hdfFile, "/");

  }

  void Close()
  {
    Flush();
    outfile_.Close();
  }


  int WriteChrSeq(const std::string& chr, const DNASequence& seq)
  {
    // Create a sequence array with the seq
    BufferedArrayPtr arr(new BufferedHDFArray<unsigned char>());
    arr->Initialize(&seq_group_.group, chr, seq.length);
    chr_to_seq_[chr] = arr;
    // Create a chromosome array
    AddChrName(chr);
    AddChrLen(seq.length);
    arr->Write((const char*)seq.seq, seq.length);

    // Add new group for track data
    HDFGroup chrgrp;
    chrgrp.Create(chr_group_, chr);
    chr_track_groups_[chr] = chrgrp;
    return 1;
  }


  // Have to AddTrack before can WriteTrackData to it
  int WriteTrackData(const TrackData& td)
  {
    std::vector<std::string> tracknames;
    track_names_.Read(tracknames);
    if (std::count(tracknames.begin(), tracknames.end(),td.track_name) == 0) {
      AddTrack(td.track_name);
    }


    return 1;
  }

  void Flush()
  {
    for (std::map<std::string, BufferedHDF2DArray<float> >::iterator it = chr_to_trackdata_.begin();
         it != chr_to_trackdata_.end(); ++it) {
      (*it).second->Flush();
    }
  }

private:
  void AddTrack(const std::string& trackname)
  {
    AddTrackName(trackname);
    std::vector<std::string> chrnames;
    chr_names_.Read(chrnames);
    // go through each chromosome group and add a new track dataset
    for (int i=0; i < chrnames.size(); ++i) {
      arr = chr_to_trackdata_[chrnames[i]];
    }

  }

  void AddTrackName(const std::string& name)
  {
    std::vector<std::string> tracknames;
    track_names_.Read(tracknames);
    tracknames.push_back(name);
    track_names_.Write(tracknames);
  }

  void AddChrName(const std::string& name)
  {
    std::vector<std::string> chrnames;
    chr_names_.Read(chrnames);
    chrnames.push_back(chrnames);
    chr_names_.Write(chrnames);
  }

  void AddChrLen(int len)
  {
    std::string s = Stringify(len);
    std::vector<std::string> chrlens;
    chr_lens_.Read(chrlens);
    chrlens.push_back(s);
    chr_lens_.Write(chrlens);
  }

  std::map<std::string, BufferedArrayPtr> chr_to_seq_;
  std::map<std::string, BufferedArrayPtr> chr_to_trackdata_;
  std::map<std::string, HDFGroup> chr_track_groups_;
  std::map<std::string, int> chr_to_len_;

  HDFFile outfile_;
  HDFAtom<std::vector<std::string> > track_names_;
  HDFAtom<std::vector<std::string> > chr_lens_;
  HDFAtom<std::vector<std::string> > chr_names_;
  std::string outfilename_;
  HDFGroup root_group_;
  HDFGroup chr_group_;
  HDFGroup seq_group_;
};

#endif
