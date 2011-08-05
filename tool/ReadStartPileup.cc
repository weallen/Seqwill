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
  BamTools::RefVector refs = b.GetReferenceData();
  TrackFile tio(trackfilename);

  LoadGenomeInfoFromChr(genome, std::string("mm9"), &g);
  std::vector<std::string> chrnames = g.chr_names();
  for (BamTools::RefVector::iterator it = refs.begin();
     it != refs.end(); ++it) {
    std::string chrname = it->RefName;
    data.plus = 0;
    data.minus = 0;
    std::cerr << "Loading chr " << chrname << std::endl;
    int refid = b.GetReferenceID(chrname);
    BamTools::BamRegion region(refid, 0, refid, g.chr_size(chrname));
    t = new Track<PlusMinusDataInt>;
    t->set_trackname(trackname);
    t->set_subtrackname(chrname);
    t->set_extends(0, g.chr_size(chrname));
    t->set_resolution(1);
    for (size_t i = 0; i < t->size(); ++i) {
      t->set(i, data);
    }
    b.SetRegion(region);
    
    int num_reads = 0;
   while (b.GetNextAlignmentCore(al)) {

     data = t->get((size_t)al.Position);
      if (al.IsReverseStrand()) {
	data.minus += 1;	
      } else {
	data.plus += 1;
      }
      t->set((size_t)al.Position, data);
    }
    b.Rewind();

    tio.WriteSubTrack<PlusMinusDataInt>(*t);
    delete t;
  }
}
