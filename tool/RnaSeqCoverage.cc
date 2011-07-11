#include <iostream>
#include <string>
#include "base/CommandLineParser.h"
#include "base/Log.h"
#include "io/BamIO.h"
#include "base/GFF.h"

int main(int argc, char** argv) {
  commandArg<std::string> iBamFilePath("-i", "BAM file path");
  commandArg<std::string> gGffPath("-g", "GFF file path");
  commandLineParser P(argc, argv);
  P.SetDescription("Compute per-exon coverage of genes in GFF file.");
  P.registerArg(iBamFilePath);
  P.registerArg(gGffPath);
  P.parse();
  std::string bamPath = P.GetStringValueFor(iBamFilePath);
  std::string gffPath = P.GetStringValueFor(gGffPath);

  std::vector<GFFelement> gff_elems;
  LoadGFFfile(gffPath, gff_elems);
  GFFelement i = gff_elems[0];
  for (std::vector<GFFelement>::iterator it = gff_elems.begin();
       it != gff_elems.end(); ++it) {
    it->Print(std::cout);
  }
  return 0;
};
