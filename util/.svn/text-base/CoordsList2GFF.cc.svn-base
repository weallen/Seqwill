//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif


#include "base/CommandLineParser.h"
#include "base/FileParser.h"


int main( int argc, char** argv )
{
  commandArg<string> cStringCmmd("-c","chromosome number e.g. IV");
  commandArg<string> lStringCmmd("-l","list of coordinates");

  commandLineParser P(argc,argv);
  P.SetDescription("Creates UCSC browser-ready GFF file from a list of coordinates");
  P.registerArg(cStringCmmd);
  P.registerArg(lStringCmmd);

  P.parse();

  string chrom = P.GetStringValueFor(cStringCmmd);
  string coords = P.GetStringValueFor(lStringCmmd);

  FlatFileParser coordsparse;
  coordsparse.Open(coords);
 
  cout << "track name=SNPs visibility=2 color=0,0,0" << endl;

  int tempcoord;

  while(coordsparse.ParseLine()) {
    tempcoord = coordsparse.AsInt(0);
    cout << "chr" << chrom << "\t.\tSNP\t" << tempcoord << "\t" << tempcoord << "\t1000\t.\t.\tchr" << chrom << endl;  
      }
   
  return 0;
}

