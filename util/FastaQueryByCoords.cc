//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif



#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{
  commandArg<string> iStringCmmd("-i","input fasta sequence");
  commandArg<string> cStringCmmd("-c","query list of coordinates to pull from fasta file");
  commandArg<bool> zStringCmmd("-zero","0-based coordinates",1);

  commandLineParser P(argc,argv);
  P.SetDescription("Queries a fasta sequence for specific coordinates and prints those bases as a list");
  P.registerArg(iStringCmmd);
  P.registerArg(cStringCmmd);
  P.registerArg(zStringCmmd);

  P.parse();

  string fasta = P.GetStringValueFor(iStringCmmd);
  string coords = P.GetStringValueFor(cStringCmmd);
  bool zero = P.GetBoolValueFor(zStringCmmd);

  FlatFileParser fastaparse;
  FlatFileParser coordsparse;
  fastaparse.Open(fasta);
  coordsparse.Open(coords);
  
  int currentPos;
  if(zero) currentPos=0;
  else currentPos=1;

  int currentCoord;

  cout << "Position   Base" << endl;

  fastaparse.ParseLine(); // consume first line of fasta file
  string templine; 
  fastaparse.GetLine(templine);

  while(coordsparse.ParseLine()) {
    currentCoord = coordsparse.AsInt(0);
      while(currentPos < currentCoord) {
	if(templine.empty()) {fastaparse.GetLine(templine);}
	templine.erase(0,1);
        currentPos++;
      }
   
  if(templine.empty()) {fastaparse.GetLine(templine);}
  cout << currentCoord << "  " << templine[0] << endl;
  templine.erase(0,1);
  currentPos++;
  
  }

  return 0;
}

