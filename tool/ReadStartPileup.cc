#include <string>
#include <vector>
#include <iostream>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

#include "base/CommandLineParser.h"
#include "base/Types.h"
#include "io/BamIO.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"

int main(int argc, char** argv) {
  commandArg<std::string> gGenome("-g", "Genome information");
  commandArg<std::string> bBam("-b", "In BAM file");
  commandArg<std::string> oDataPath("-o", "Trackfile name");
  commandArg<std::string> tTrackName("-t", "Out trackname");
  commandLineParser P(argc, argv);
  P.SetDescription("Save pileup of read starts from BAM");
  P.registerArg(gGenome);
  P.registerArg(bBam);
  P.registerArg(oDataPath);
  P.registerArg(tTrackName);
  P.parse();

  std::string genome = P.GetStringValueFor(gGenome);
  std::string bamname = P.GetStringValueFor(bBam);
  std::string trackfilename = P.GetStringValueFor(oDataPath);
  std::string trackname = P.GetStringValueFor(tTrackName);

  Track<PlusMinusDataInt>* t;
  GenomeInfo g;
  BamTools::BamAlignment al;
  BamTools::BamReader b;
  PlusMinusDataInt data;
  data.plus = 0;
  data.minus = 0;
  if (!b.Open(bamname)) {
    std::cerr << "Couldn't open input BAM file." << std::endl;
    return -1;
  }
  TrackFile tio(trackfilename);

  LoadGenomeInfoFromChr(genome, std::string("mm9"), &g);
  std::vector<std::string> chrnames = g.chr_names();
  for (std::vector<std::string>::iterator it = chrnames.begin();
       it != chrnames.end(); ++it) {
    std::cerr << "Loading chr " << *it << std::endl;
    int refid = b.GetReferenceID(*it);
    BamTools::BamRegion region(refid, 0, refid, g.chr_size(*it));
    t = new Track<PlusMinusDataInt>;
    t->set_trackname(trackname);
    t->set_subtrackname(*it);
    t->set_extends(0, g.chr_size(*it));
    t->set_resolution(1);
    for (size_t i = 0; i < t->size(); ++i) {
      t->set(i, data);
    }
    b.SetRegion(region);
    while (b.GetNextAlignmentCore(al)) {
      data = t->get((int)al.Position);
      if (al.IsReverseStrand()) {
	data.minus += 1;
      } else {
	data.plus += 1;
      }
      t->set((int)al.Position, data);
    }
    tio.WriteSubTrack<PlusMinusDataInt>(*t);
    delete t;
  }
}
