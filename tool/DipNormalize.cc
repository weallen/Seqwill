#include <string>
#include <vector>
#include <iostream>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

#include <omp.h>

#include "base/CommandLineParser.h"
#include "base/FileUtil.h"
#include "base/Types.h"
#include "io/BamIO.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"
#include "analysis/MedipNormalize.h"

int main(int argc, char** argv) {
   
  omp_set_num_threads(25);

  commandArg<std::string> gGenome("-g", "Genome information");
  commandArg<std::string> bBam("-b", "In BAM file");
  commandArg<std::string> cChr("-c", "Chromosome name");
  commandArg<std::string> oDataPath("-o", "Trackfile name");
  commandArg<std::string> tTrackName("-t", "Out trackname");
  commandArg<int> nRes("-n", "Output resolution", 100);
  
  commandArg<int> eExtend("-e", "Extend length if not paired end", -1);
  

  commandLineParser P(argc, argv);
  P.SetDescription("Normalize BAM of me or hmeDIP reads by CpG content.");
  P.registerArg(gGenome);
  P.registerArg(bBam);
  P.registerArg(oDataPath);
  P.registerArg(tTrackName);
  P.registerArg(cChr);
  P.registerArg(nRes);
  P.registerArg(eExtend);
  P.parse();

  std::string genome = P.GetStringValueFor(gGenome);
  std::string bamname = P.GetStringValueFor(bBam);
  std::string trackfilename = P.GetStringValueFor(oDataPath);
  std::string trackname = P.GetStringValueFor(tTrackName);
  std::string chrname = P.GetStringValueFor(cChr);
  int res = P.GetIntValueFor(nRes);
  int extend = P.GetIntValueFor(eExtend);

  Track<float> t;
  GenomeInfo g;

  if (!FileExists(bamname)) {
    std::cerr << "Couldn't find bam file " << bamname << std::endl;
    return -1;
  }
   
  if (!FileExists(genome)) {
    std::cerr << "Couldn't find genome file " << genome << std::endl;
    return -1;
  }

  BamIO* b = new BamIO(bamname);
  SingleReadFactory* reads = b->LoadChrSingleReads(chrname);

  BamTools::RefVector refs = b->reader()->GetReferenceData();
  bool has_chr = false;
  for (BamTools::RefVector::iterator it = refs.begin();
       it != refs.end(); ++it) {
    if (it->RefName == chrname) 
      has_chr = true;
  }
  if (!has_chr) {
    std::cerr << "Couldn't find " << chrname << " in bam data" << std::endl;
    return -1;
  }

  // Load chr
  LoadGenomeInfoFromChr(genome, std::string("mm9"), &g);
  std::vector<std::string> chrnames = g.chr_names();
  if (std::find(chrnames.begin(), chrnames.end(), chrname) == chrnames.end()) {
    std::cerr << "Couldn't find chr " << chrname << " in reference sequences" << std::endl;
    return -1;
  }
  TrackFile chrio(genome);
  Track<unsigned char>::Ptr chr(new Track<unsigned char>);
  chrio.ReadSubTrack<unsigned char>(std::string("mm9"), chrname, *chr);
  
  // finally actually do the computation
  TrackFile tio(trackfilename);
  
  MedipNormalize norm;
  if (extend > -1) {
    norm.set_frag_len(extend);
  }
  norm.set_chr(chr);
  norm.set_reads(reads);
  norm.set_resolution(res);
  norm.set_out_track_name(trackname);
  norm.set_out_subtrack_name(chrname);
  norm.Compute();
  Track<float>::Ptr out = norm.output();
  tio.WriteSubTrack<float>(*out);

  delete reads;
  delete b;
  return 1;
}
