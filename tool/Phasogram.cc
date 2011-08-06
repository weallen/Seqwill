#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

#include "base/CommandLineParser.h"
#include "base/Types.h"
#include "base/Interval.h"
#include "base/BEDelement.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"
#include "analysis/Nucleosome.h"

int main(int argc, char** argv) {
  commandArg<std::string> iTrackfile("-i", "In track file");
  commandArg<std::string> tTrackName("-t", "In trackname");
  commandArg<std::string> oDataPath("-o", "Phasogram filename");
  commandArg<int> nThreshold("-n", "Threshold");
  commandArg<int> lLength("-l", "Histogram size");
  commandArg<std::string> bBed("-b", "BED file of regions to compute histogram over");

  commandLineParser P(argc, argv);
  P.SetDescription("Save pileup of read starts from BAM");
  P.registerArg(iTrackfile);
  P.registerArg(tTrackName);
  P.registerArg(oDataPath);
  P.registerArg(nThreshold);
  P.registerArg(lLength);
  P.registerArg(bBed);
  P.parse();

  std::string datapath = P.GetStringValueFor(oDataPath);
  std::string trackfilename = P.GetStringValueFor(iTrackfile);
  std::string trackname = P.GetStringValueFor(tTrackName);
  int threshold = P.GetIntValueFor(nThreshold);
  int length = P.GetIntValueFor(lLength);
  std::string bedfile = P.GetStringValueFor(bBed);

  Track<PlusMinusDataInt>* t;
  TrackFile tio(trackfilename);
  PlusMinusDataInt data;
  Histogram h(length);
  std::map<std::string, std::vector<BEDelement> > chrs;
  std::vector<BEDelement> bedelems;
  std::cerr << "Loading bedfile " << bedfile << std::endl;
  LoadBEDfile(bedfile, bedelems);
  
  for (std::vector<BEDelement>::iterator it = bedelems.begin();
       it != bedelems.end(); ++it) {
    chrs[it->GetName()].push_back(*it);
  }
  std::cerr << "Computing phasogram" << std::endl;
  std::vector<std::string> chrnames = tio.GetSubTrackNames(trackname);  
  for (std::vector<std::string>::iterator it = chrnames.begin();
       it != chrnames.end(); ++it) {
    std::string chrname = *it;
    std::cerr << chrname << std::endl;
    bedelems = chrs[chrname];
    t = new Track<PlusMinusDataInt>();
    tio.ReadSubTrack<PlusMinusDataInt>(trackname, chrname, *t);
    for (std::vector<BEDelement>::iterator it = bedelems.begin();
	 it != bedelems.end(); ++it) {      
      for (int i =it->GetStart() ; i < (it->GetEnd() - length); ++i) {
	data = t->get(i);
	if (data.plus > threshold) {
	  for (int j = 0; j < length; ++j) {
	    data = t->get(i + j);
	    h.add_to_bin(j, data.plus);
	  }
	}
      }
      delete t;
    }
  }
  std::fstream outfile;
  outfile.open(datapath.c_str(), std::fstream::out);
  for (int i = 0; i < h.num_bins(); ++i) {
    outfile << h.get_bin(i) << std::endl;
  }    
}
