#include <iostream>
#include <string>
#include <H5Cpp.h>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "analysis/Genome.h"

using namespace std;

int main(int argc, char** argv) {
  commandArg<string> iGenomeFasta("-i", "genome fasta dir");
  commandArg<string> aMedip("-a", "first medip bam");
  commandArg<string> bMedip("-b", "second medip bam");
  
  commandLineParser P(argc, argv);
  P.SetDescription("Comparse two medips bam files with an HMM.");
  P.registerArg(iGenomeFasta);
  P.registerArg(aMedip);
  P.registerArg(bMedip);

  P.parse();

  string genomeFasta = P.GetStringValueFor(iGenomeFasta);
  string medip1 = P.GetStringValueFor(aMedip);
  string medip2 = P.GetStringValueFor(bMedip);


  cout << "Reading genome fasta directory " << genomeFasta << "..." << endl;
  
  // Read genome and build scaffold
  cout << "Reading first medip " << medip1 << "..." << endl;

  cout << "Reading second medip " << medip2 << "..." << endl;
  return 0;
}

