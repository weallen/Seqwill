#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/StringUtil.h"

#include "analysis/Chromosome.h"

using namespace std;

int main(int argc, char** argv) {
  commandArg<string> iGenomeFasta("-i", "genome fasta dir");
  commandArg<string> oDataPath("-o", "HDF5 dir name"); 
  commandArg<string> gGenome("-g", "Genome name");
  commandLineParser P(argc, argv);
  P.SetDescription("Load fasta sequence files into HDF5.");
  P.registerArg(iGenomeFasta);
  P.registerArg(oDataPath);
  P.registerArg(gGenome);
  P.parse();

  string genomeFasta = P.GetStringValueFor(iGenomeFasta);
  string dataPath = P.GetStringValueFor(oDataPath);
  string genomeName = P.GetStringValueFor(gGenome);
  cout << "Reading genome fasta directory " << genomeFasta << "..." << endl;

  DIR* fastaDir;
  struct dirent* dirp; 

  fastaDir = opendir(genomeFasta.c_str());
  if (fastaDir == NULL) {
    cerr << "Fasta dir " << genomeFasta << " not found." << endl;
    return -1;
  }

  while ((dirp = readdir(fastaDir))) {
    string fname(dirp->d_name);
    string extension = ".fa";
    if (Contains(fname, extension)) {
      cout << "Processing " << fname << endl; 
      SaveChrFromFASTA(dataPath, genomeFasta + "/" + fname, genomeName);
    }
  }

  closedir(fastaDir);
  return 0;
}

