#include <string>
#include <vector>
#include <iostream>


#include "base/CommandLineParser.h"
#include "base/FileUtil.h"
#include "base/Types.h"
#include "io/TrackIO.h"

int main(int argc, char** argv) {
   
  commandArg<std::string> iDataPath("-i", "In trackfile");
  commandArg<std::string> oDataPath("-o", "Trackfile name"); 

  commandLineParser P(argc, argv);
  P.SetDescription("Normalize track of counts by RPM");

  P.registerArg(oDataPath);
  P.registerArg(iDataPath);
  P.parse();

  std::string outtrackfilename = P.GetStringValueFor(oDataPath);
	std::string intrackfilename = P.GetStringValueFor(iDataPath);


  if (!FileExists(intrackfilename)) {
    std::cerr << "Couldn't find track file " << intrackfilename << std::endl;
    return -1;
  }
   

  TrackFile tin(intrackfilename);
	TrackFile tout(outtrackfilename);

  std::vector<std::string> curr_tracknames = tin.GetTrackNames();
  for (std::vector<std::string>::const_iterator tname = curr_tracknames.begin();
      tname != curr_tracknames.end(); ++tname) {
			std::cout << "Processing track " << *tname << std::endl;
      std::vector<std::string> curr_subtracknames = tin.GetSubTrackNames(*tname);
			float sum = 0.0;
      for (std::vector<std::string>::const_iterator stname = curr_subtracknames.begin();
          stname != curr_subtracknames.end(); ++stname) { 
          Track<float>::Ptr track(new Track<float>);
          tin.ReadSubTrack<float>(*tname, *stname, *track);			
					for (size_t i = 0; i < track->size(); ++i) {
							sum += track->get(i);
					}
      }
      for (std::vector<std::string>::const_iterator stname = curr_subtracknames.begin();
          stname != curr_subtracknames.end(); ++stname) { 
          Track<float>::Ptr track(new Track<float>);
          tin.ReadSubTrack<float>(*tname, *stname, *track);
					for (size_t i = 0; i < track->size(); ++i) {
							float temp = (1.0e6 * track->get(i)) / sum;
							track->set(i, temp);
					}	
					tout.WriteSubTrack<float>(*track);
      } 
  }

  
  return 1;
}
