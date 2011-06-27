#ifndef HDFTRACKREADER_H_
#define HDFTRACKREADER_H_

#include <iostream>
#include <string>
#include <vector>

#include "hdf5/BufferedHDF2DArray.h"
#include "hdf5/HDFGroup.h"
#include "hdf5/HDFAtom.h"
#include "hdf5/HDFFile.h"

#include "data/TrackData.h"
#include "data/ChrData.h"
#include "data/SeqData.h"

/* Use BufferedHDFArray to read
 * contiguous blocks of data from each track
 * in parallel.
 *
 * Load a chromosome at a time.
 *
 */


class HDFTrackReader
{
public:
  HDFTrackReader()
    : curr_chr_(0)
    , ntracks_(0)
    {}

  virtual ~HDFTrackReader() {}

  int Init(const std::vector<std::string>& genomefilename) {
    try {
        genome_file_.Open(genomefilename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    } catch (Exception& e) {
      std::cout << e.getDetailMsg() << std::endl;
      return 0;
    }
    if (track_names_.Initialize(genome_file_, 'TrackNames') == 0) {
      return 0;
    }
    if (chr_names_.Initialize(genome_file_, 'ChrNames') == 0) {
      return 0;
    }
  }

  void Open() {
    try {
        genome_file_.Open(genomefilename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    } catch (Exception& e) {
      std::cout << e.getDetailMsg() << std::endl;
      return 0;
    }
  }

  void Close() {
    genome_file.close();
  }


  void ReadChrData(const std::string& chrname, ChrData* chrdata)
  {
    track_names_.Read(chrdata->tracknames);
    chr_names_.Read(chrdata->chrnames);
    chr_lens_.Read(chrdata->chrlens);
  }

  void ReadChrSeq(const std::string& chrname, SeqData* sd)
  {
    sd->chr_name = chrname;
    ChrData cr;
    GetChrData(chrname, &cr);

  }

  void ReadChrTrack(const std::string& chrname, const std::string& trackname, TrackData* td)
  {
  }

  void ReadChrRegionTrack(const std::string& trackname, const std::string& chrname,
                         int start, int end, TrackData* td)
  {
  }
private:
  HDFFile genome_file_;
  std::vector<HDFGroup> chr_groups_;
  BufferedHDF2DArray<float> trackdata_;
  HDFAtom<std::vector<std::string> > track_names_;
  HDFAtom<std::vector<std::string> > chr_names_;
  HDFAtom<std::vector<int> > chr_lens_;
  int ntracks_;
  int curr_chr_;

};
#endif
