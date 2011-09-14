#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "base/CommandLineParser.h"
#include "base/Types.h"
#include "base/BEDelement.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"
#include "analysis/MedipNormalize.h"


int main(int argc, char** argv) {
  commandArg<std::string> pPosition("-i", "Nucleosome position trackfile");
  commandArg<std::string> tTrack("-t", "Nucleosome position track name");
  commandArg<std::string> bBed("-b", "BED file of regions");
  commandArg<int> fFwd("-f","Distance to extend forward from each region start",1000);
  commandArg<int> rRev("-r","Distance to extend backward from each region start",1000);
  commandArg<std::string> oDataPath("-o", "Out datafile name");

  commandLineParser P(argc, argv);
  P.SetDescription("Get nucleosome positioning for regions");
  P.registerArg(pPosition);
  P.registerArg(tTrack);
  P.registerArg(bBed);
  P.registerArg(fFwd);
  P.registerArg(rRev);
  P.registerArg(oDataPath);
  P.parse();

  std::string trackfilename = P.GetStringValueFor(pPosition);
  std::string trackname = P.GetStringValueFor(tTrack);
  std::string bedfile = P.GetStringValueFor(bBed);
  int fwd = P.GetIntValueFor(fFwd);
  int rev = P.GetIntValueFor(rRev);
  std::string outname = P.GetStringValueFor(oDataPath);

  Track<float>* t;

  TrackFile tio(trackfilename);

  std::map<std::string, std::vector<BEDelement> > chrs;
  std::vector<BEDelement> bedelems;

  std::cerr << "Loading bedfile " << bedfile << std::endl;
  LoadBEDfile(bedfile, bedelems);  

  for (std::vector<BEDelement>::iterator it = bedelems.begin();
       it != bedelems.end(); ++it) {
    chrs[it->GetChromosome()].push_back(*it);
  }

  std::fstream outfile;
  outfile.open(outname.c_str(), std::fstream::out);

  std::cerr << "Computing positioning" << std::endl;
  std::vector<std::string> chrnames = tio.GetSubTrackNames(trackname);  

  for (std::vector<std::string>::iterator it = chrnames.begin();
       it != chrnames.end(); ++it) {
    std::string chrname = *it;
    std::cerr << chrname << std::endl;
    bedelems = chrs[chrname];
    t = new Track<float>;
    tio.ReadSubTrack<float>(trackname, chrname, *t);
    for (size_t i = 0; i < bedelems.size(); ++i) {
      BEDelement b = bedelems[i];
      int start = b.GetStart() - rev;
      int stop = b.GetStart() + fwd;
      for (int j = start; j < stop; ++j) {
	outfile << (t->get(j)) << " ";
      }
      outfile << std::endl;
    }
    delete t;
  }
  outfile.close();
}
