#include <iostream>
#include <string>
#include <sys/dir.h>
#include <sys/types.h>
#include <sys/param.h>
#include <stdio.h>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "analysis/Genome.h"

using namespace std;

int main(int argc, char** argv) {
  commandArg<string> iDataPath("-i", "HDF5 dir");
  commandArg<string> tTrack("-t", "data bam file");
  commandLineParser P(argc, argv);
  P.SetDescription("Load sequence data into HDF5.");
  P.registerArg(iDataPath);
  P.registerArg(tTrack);
  P.parse();

  string dataPath = P.GetStringValueFor(iDataPath);
  string trackPath = P.GetStringValueFor(tTrack);

  TrackReader g(dataPath);

  cout << "Reading bam file" << trackPath << "..." << endl;
  
  return 0;
}

