#ifndef HDFTRACKREADER_H_
#define HDFTRACKREADER_H_

#include <iostream>
#include <string>
#include <vector>

#include "hdf5/BufferedHDF2DArray.h"
#include "hdf5/HDFGroup.h"
#include "hdf5/HDFAtom.h"
#include "hdf5/HDFFile.h"

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
  HDFTrackReader() {}
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

  void Close() {
    genome_file.close();
  }

private:
  HDFFile genome_file_;
  HDFChrReader chr_reader_;
  HDFAtom<std::vector<std::string> > track_names_;
  HDFAtom<std::vector<std::string> > chr_names_;
  int ntracks_;

};
#endif
