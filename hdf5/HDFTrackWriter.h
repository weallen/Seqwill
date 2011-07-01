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


// TRACKFILE :
//
//
template <typename TypeT>
class HDFTrackWriter
{
public:

  void Create(const char* outfilename) {
    std::string s(outfilename);
    Create(s);
  }

  void Create(std::string& outfilename) {
    std::vector<std::string> empty_string;
    std::string s;

    outfilename_ = outfilename;
    outfile_.Create(outfilename_);
    root_group_.Initialize(*outfile_.hdfFile, "/");
    metadata_.Create(root_group_.group, "Metadata");
    metadata_.Write(s);
    subtrack_names_.Create(root_group_.group, "Subtracks");
    subtrack_names_.Write(empty_string);
  }

  int Open(const char* outfilename)
  {
    std::string s(outfilename);
    return Open(s);
  }

  int Open(std::string& outfilename)
  {
    std::vector<std::string> subtracknames;
    outfilename_ = outfilename;
    try {
      outfile_.Create(outfilename);
    } catch (Exception& e) {
      std::cout << e.getDetailMsg() << std::endl;
      return -1;
    }
    root_group_.Initialize(*outfile_.hdfFile, "/");
    subtrack_names_.Initialize(root_group_.group, "Subtracks");
    metadata_.Initialize(root_group_.group, "Metadata");
    subtrack_names_.Read(subtracknames);
    for (std::vector<std::string>::iterator i = subtracknames.begin();
         i != subtracknames.end(); ++i) {
      BufferedHDFArray<TypeT>* currtrack = new BufferedHDFArray<TypeT>();
      currtrack->Initialize(&root_group_.group, *i);
      name_to_subtracks_[*i] = currtrack;
    }
    return 1;
  }

  void Close()
  {
    Flush();
    outfile_.Close();
    root_group_.Close();
    subtrack_names_.Close();
    metadata_.Close();
  }

  int WriteSubTrackData(int start, int end, const TrackData<TypeT>& td)
  {
    BufferedHDFArray<TypeT>* arr;
    std::vector<TypeT> dset;
    int tdlen = end - start;

    if (name_to_subtracks_.count(td.track_name) == 0) {
      return WriteSubTrackData(td);
    }
    if ((td.end - td.start) != (end - start)) {
      return -1;
    }
    arr = name_to_subtracks_[td.track_name];
    arr->ReadDataset(dset);
    std::copy(td.track, td.track + tdlen, dset[td.start]);
    arr->Write(dset);
    return 1;
  }

  // Have to AddTrack before can WriteTrackData to it
  int WriteSubTrackData(const TrackData<TypeT>& td)
  {
    std::vector<std::string> tracknames;
    BufferedHDFArray<TypeT>* arr = new BufferedHDFArray<TypeT>();
    HDFAtom<int> subtracklen_atom;
    int tdlen = td.end - td.start;
    int subtracklen;
    HDFAtom<std::string> currtrackname_atom;

    subtrack_names_.Read(tracknames);
    if (std::count(tracknames.begin(), tracknames.end(),td.track_name) == 0) {
      AddSubTrackName(td.track_name);
      arr->Initialize(&root_group_.group, td.track_name, td.end - td.start);
      subtracklen_atom.Create(*arr, "Length");
      subtracklen_atom.Write(td.end - td.start);
      currtrackname_atom.Initialize(*arr, "Name");
      currtrackname_atom.Write(td.track_name);
      name_to_subtracks_[td.track_name] = arr;
    } else {
      arr = name_to_subtracks_[td.track_name];
      subtracklen_atom.Initialize(*arr, "Length");
      subtracklen_atom.Read(subtracklen);
      if (subtracklen != tdlen) {
        return -1;
      } else {
        arr->Write(td.track, tdlen);
      }
    }
    return 1;
  }

  void Flush()
  {
    std::vector<std::string> sn = GetSubtrackNames();
    for (std::vector<std::string>::iterator it = sn.begin();
         it != sn.end(); ++it) {
        name_to_subtracks_[*it]->Flush();
    }
  }

  // Store metadata in protobuf serialized format
  void AddMetadata(const std::string& data)
  {
    metadata_.Write(data);
  }

  const std::string GetMetadata()
  {
    std::string mdata;

    metadata_.Read(mdata);
    return mdata;
  }

  const std::vector<std::string> GetSubtrackNames()
  {
    std::vector<std::string> subtracknames;
    subtrack_names_.Read(subtracknames);
    return subtracknames;
  }

private:

  void AddSubTrackName(const std::string& name)
  {
    std::vector<std::string> tracknames;
    subtrack_names_.Read(tracknames);
    tracknames.push_back(name);
    subtrack_names_.Write(tracknames);
  }


  HDFFile outfile_;
  std::string outfilename_;

  HDFGroup root_group_;
  HDFAtom<std::string> metadata_;
  std::map<std::string, BufferedHDFArray<TypeT>* > name_to_subtracks_;
  HDFAtom<std::vector<std::string> > subtrack_names_;
};

#endif
