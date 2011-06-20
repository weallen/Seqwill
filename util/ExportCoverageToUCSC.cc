
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

int main(int argc, char** argv)
{

  commandArg<string> covFile("-f","file of coverages");
  commandArg<string> chr("-c","target chromosome");
  commandArg<int> targId("-t","target id");
  commandArg<int> offset("-o","genomic coord for start of target");
  commandArg<string> trackname("-n","name for track in UCSC");

  commandLineParser P(argc,argv);

  P.registerArg(covFile);
  P.registerArg(targId);
  P.registerArg(offset);
  P.registerArg(chr);
  P.registerArg(trackname);
  
  P.parse();

  string some = P.GetStringValueFor(chr);
  string file  = P.GetStringValueFor(covFile);
  int id = P.GetIntValueFor(targId);
  int start = P.GetIntValueFor(offset);
  string name(P.GetStringValueFor(trackname));

  FlatFileParser fp(file);
  string line("");

//   cout << "track type=wiggle_0" << endl;
//   cout << "variableStep chrom=" << some << endl;

  cout << "track type=bedGraph name=" << flush; 
  cout << name << endl;
  int count(0);
  while ( fp.GetLine(line) )
  {
    StringParser sp;
    sp.SetLine(line);

    if ( sp.AsInt(0) != id )
      continue;

    cout << some <<"\t" << sp.AsInt(1)+start 
	 <<"\t"
	 << sp.AsInt(2)+start 
	 <<"\t"
	 << sp.AsInt(3)
	 << endl;

  }
}
