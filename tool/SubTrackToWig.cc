#include <vector>
#include <fstream>
#include <iostream>
#include "common/Track.h"
#include "io/TrackIO.h"
#include "base/CommandLineParser.h"

// Note: for now just writes PlusMinusData
int main(int argc, char** argv) {
  commandArg<std::string> iTrackfile("-i", "In track file");
  commandArg<std::string> tTrackName("-t", "Trackname");
  commandArg<std::string> sSubTrack("-s", "Subtrackname");
  commandArg<std::string> oWig("-o", "WIG filename");

  commandLineParser P(argc, argv);
  P.SetDescription("Save track data as wig");
  P.registerArg(iTrackfile);
  P.registerArg(tTrackName);
  P.registerArg(sSubTrack);
  P.registerArg(oWig);
  P.parse();

  std::string trackfile = P.GetStringValueFor(iTrackfile);
  std::string track = P.GetStringValueFor(tTrackName);
  std::string subtrack = P.GetStringValueFor(sSubTrack);
  std::string wig = P.GetStringValueFor(oWig); 

  Track<float> t;
  TrackFile tio(trackfile);

  std::vector<std::string> tnames = tio.GetTrackNames();
  if (std::find(tnames.begin(), tnames.end(), track) == tnames.end()) {
    std::cerr << "Couldn't find track " << track << std::endl;
    return -1;
  }

  std::vector<std::string> stnames = tio.GetSubTrackNames(track);
  if (std::find(stnames.begin(), stnames.end(), subtrack) == stnames.end()) {
    std::cerr << "Couldn't find subtrack " << subtrack << std::endl;
    return -1;
  }

  tio.ReadSubTrack<float>(track, subtrack, t);
  float d;

  std::fstream plus;
  
  plus.open((wig + ".wig").c_str(), std::fstream::out);
//  minus.open((wig + "_minus.wig").c_str(), std::fstream::out);

  plus << "fixedStep chrom=" << subtrack  << " start=0 step=" << t.resolution() << std::endl;
//  minus << "fixedStep chrom=" << subtrack  << " start=0 step=" << t.resolution() << std::endl;
  for (size_t i = 0; i < t.size(); ++i) {
    d = t.get(i);
    plus << d << std::endl;
  //  minus << d.minus << std::endl;
  }
  plus.close();
//  minus.close();
  return 1;
}
