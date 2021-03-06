#include <string>
#include <vector>
#include <iostream>

#include "base/CommandLineParser.h"
#include "base/Types.h"
#include "io/BamIO.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"
#include "analysis/MedipNormalize.h"


int main(int argc, char** argv) {
  commandArg<std::string> gGenome("-g", "Genome information");
  commandArg<std::string> oDataPath("-o", "Out trackfile name");
  commandArg<int> nRes("-n", "Resolution", 100);

  commandLineParser P(argc, argv);
  P.SetDescription("Count CpGs in genome.");
  P.registerArg(gGenome);
  P.registerArg(oDataPath);
  P.registerArg(nRes);
  P.parse();

  std::string genome = P.GetStringValueFor(gGenome);
  std::string trackfilename = P.GetStringValueFor(oDataPath);
  int res = P.GetIntValueFor(nRes);

  Track<int> t;
  GenomeInfo g;

  TrackFile chrio(genome);
  TrackFile tio(trackfilename);

  LoadGenomeInfoFromChr(genome, std::string("mm9"), &g);
  std::vector<std::string> chrnames = g.chr_names();

	CpGCounter* counter = new CpGCounter;

	counter->set_out_track_name(std::string("mm9"));
	counter->set_resolution(res);
  for (std::vector<std::string>::iterator it = chrnames.begin();
       it != chrnames.end(); ++it) {
    Track<unsigned char>::Ptr chr(new Track<unsigned char>);
    if (it->find(std::string("random")) == std::string::npos) {
	    std::cerr << "Processing chromosome " << *it << std::endl;
	    chrio.ReadSubTrack<unsigned char>(std::string("mm9"), *it, *chr);
	    counter->set_input(chr);
	    counter->set_out_subtrack_name(*it);
	    counter->Compute();
	    Track<int>::Ptr out = counter->output();
	    tio.WriteSubTrack<int>(*out);
	    //chr.reset();
    }
  }

  delete counter;
}
