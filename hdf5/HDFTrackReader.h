#ifndef HDFTRACKREADER_H_
#define HDFTRACKREADER_H_

#include <iostream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "hdf5/BufferedHDFArray.h"
#include "hdf5/HDFGroup.h"
#include "hdf5/HDFAtom.h"
#include "hdf5/HDFFile.h"

#include "data/TrackData.h"
#include "data/ChrData.h"
#include "data/SeqData.h"



template <typename TypeT>
class HDFTrackReader
{
public:
  HDFTrackReader() {}

  virtual ~HDFTrackReader() {}

  int Open(const std::string& filename) {
    filename_ = filename;
    std::vector<std::string> subtracknames;

    try {
        file_.Create(filename);
    } catch (Exception& e) {
      std::cout << e.getDetailMsg() << std::endl;
      return -1;
    }
    root_group_.Initialize(*file_.hdfFile, "/");
    if (subtrack_names_.Initialize(root_group_, "Subtracks") == 0) {
      return -1;
    }
    subtrack_names_.Read(subtracknames);
    for (std::vector<std::string>::iterator i = subtracknames.begin();
         i != subtracknames.end(); ++i) {
      BufferedHDFArray<TypeT>* arr = new BufferedHDFArray<TypeT>();
      arr->Initialize(&root_group_.group, *i);
      name_to_subtracks_[*i] = arr;
    }
  }

  void Close() {
    file_.Close();
    root_group_.Close();
    std::vector<std::string> sn;
    subtrack_names_.Read(sn);
    for (std::vector<std::string>::iterator i = sn.begin(); i != sn.end(); ++i) {
      name_to_subtracks_[*i]->Close();
      delete name_to_subtracks_[*i];
    }
    subtrack_names_.Close();
    metadata_.Close();
  }

  std::string GetMetadata()
  {
    std::string mdata;
    metadata_.Read(mdata);
    return mdata;
  }

  std::vector<std::string> GetSubtrackNames()
  {
    std::vector<std::string> subtracknames;
    subtrack_names_.Read(subtracknames);
    return subtracknames;
  }

  std::vector<int> GetSubtrackLengths()
  {
    std::vector<int> subtracklens;
    std::vector<std::string> subtracknames = GetSubtrackNames();
    HDFAtom<int> currlen_atom;
    int currlen;
    BufferedHDFArray<TypeT>* currarr;

    for (std::vector<std::string>::iterator i = subtracknames.begin();
         i != subtracknames.end(); ++i) {
      currarr = name_to_subtracks_[*i];
      currlen_atom.Initialize(*currarr, "Length");
      currlen_atom.Read(currlen);
      subtracklens.push_back(currlen);
    }
    return subtracklens;
  }

  int ReadSubtrackRegion(const std::string& trackname, int start, int end, TrackData<TypeT>* td)
  {
    int tlen = end - start;
    HDFAtom<int> subtracklen_atom;
    int subtracklen;
    BufferedHDFArray<TypeT>* arr;

    if (name_to_subtracks_.count(trackname) == 0) {
      return -1;
    }
    arr = name_to_subtracks_[trackname];
    subtracklen_atom.Initialize(*arr, "Length");
    subtracklen_atom.Read(subtracklen);
    if (end > subtracklen) {
      return -1;
    }
    td->track = new TypeT[tlen];
    arr->Read(start, end, td->track);
    td->track_name = trackname;
    td->start = start;
    td->end = end;
  }

  int ReadSubtrack(const std::string& trackname, TrackData<TypeT>* td)
  {
    BufferedHDFArray<TypeT>* arr;
    HDFAtom<int> subtracklen_atom;
    int subtracklen;

    if (name_to_subtracks_.count(trackname) == 0) {
      return -1;
    }
    arr = name_to_subtracks_[trackname];
    subtracklen_atom.Initialize(*arr, "Length");
    subtracklen_atom.Read(subtracklen);

    return ReadSubtrackRegion(trackname, 0, subtracklen, td);
  }

private:

  HDFFile file_;
  std::string filename_;

  HDFGroup root_group_;
  HDFAtom<std::string> metadata_;
  std::map<std::string, BufferedHDFArray<TypeT>* > name_to_subtracks_;
  HDFAtom<std::vector<std::string> > subtrack_names_;
};
#endif
