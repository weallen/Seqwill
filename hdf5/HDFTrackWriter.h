#ifndef HDFTRACKWRITER_H_
#define HDFTRACKWRITER_H_

#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include "hdf5/HDFGroup.h"
#include "hdf5/HDFFile.h"
#include "hdf5/HDFAtom.h"
#include "hdf5/BufferedHDF2DArray.h"
#include "hdf5/BufferedHDFArray.h"

#include "base/StringUtil.h"
#include "base/DNASequence.h"

#include "data/TrackData.h"

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

  void Open(std::string& outfilename) {
    outfilename_ = outfilename;
    outfile_.Create(outfilename_);
    root_group_.Initialize(*outfile_.hdfFile, "/");

  }

  void Close() {
    Flush();
    outfile_.Close();
  }


  int WriteChrSeq(const std::string& chr, const DNASequence& seq) {

    // Create a sequence array with the seq
    BufferedHDFArray<unsigned char> seq_array;
    seq_array.Initialize(&seq_group_.group, chr, seq.length);
    chr_to_seq_[chr] = seq_array;
    // Create a chromosome array
    AddChrName(chr);
    AddChrLen(seq.length);
    return 1;
  }

  void AddTrack(const std::string& trackname) {
    AddTrackName(trackname);
    std::vector<std::string> chrnames;
    chr_names_.Read(chrnames);
    // go through each chromosome file and add a row
    for (std::vector<std::string>::const_iterator it = chrnames.begin();
         it != chrnames.end(); ++it) {
   //   chr_to_trackdata_[*it].
    }
  }

  // Have to AddTrack before can WriteTrackData to it
  int WriteTrackData(const TrackData& td) {
    std::vector<std::string> tracknames;
    track_names_.Read(tracknames);
    if (std::count(tracknames.begin(), tracknames.end(),td.track_name) == 0) {
      return 0;
    }

    return 1;
  }

  void Flush() {
    for (std::map<std::string, BufferedHDF2DArray<float> >::iterator it = chr_to_trackdata_.begin();
         it != chr_to_trackdata_.end(); ++it) {
      (*it).second->Flush();
    }
  }

  private:
  void AddTrackName(const std::string& name) {
    std::vector<std::string> tracknames;
    track_names_.Read(tracknames);
    tracknames.push_back(name);
    track_names_.Write(tracknames);
  }

  void AddChrName(const std::string& name) {
    std::vector<std::string> chrnames;
    chr_names_.Read(chrnames);
    chrnames.push_back(chrnames);
    chr_names_.Write(chrnames);
  }

  void AddChrLen(int len) {
    std::string s = Stringify(len);
    std::vector<std::string> chrlens;
    chr_lens_.Read(chrlens);
    chrlens.push_back(s);
    chr_lens_.Write(chrlens);
  }
  std::map<std::string, BufferedHDFArray<unsigned char> > chr_to_seq_;
  std::map<std::string, BufferedHDF2DArray<float> > chr_to_trackdata_;
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
