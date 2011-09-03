#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <vector>



#include "analysis/Genome.h"
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/StringUtil.h"
#include "analysis/Nucleosome.h"
#include "io/Traits.h"
#include "common/Track.h"
#include "io/TrackIO.h"

using namespace std;

int main(int argc, char** argv) {
    commandArg<string> tTrackfile("-i", "In trackfile");
    commandArg<string> oTrackfile("-o", "Out trackfile");
    commandArg<string> tTrack("-t", "Track name");

    commandLineParser P(argc, argv);
    P.SetDescription("Find nucleosome positioning stringency.");
    P.registerArg(tTrackfile);
    P.registerArg(oTrackfile);
    P.registerArg(tTrack);
    P.parse();
    

    string trackname = P.GetStringValueFor(tTrack);
    string in_trackfile = P.GetStringValueFor(tTrackfile);
    string out_trackfile = P.GetStringValueFor(oTrackfile);
    
    
    TrackFile in_file(in_trackfile);
    TrackFile out_file(out_trackfile);
    std::vector<std::string> tnames = in_file.GetTrackNames();
    if (std::find(tnames.begin(), tnames.end(), trackname) == tnames.end()) {
        std::cerr << "Couldn't find track " << trackname << std::endl;
        return -1;
    }

    std::vector<std::string> curr_subtracknames = in_file.GetSubTrackNames(trackname);
    for (std::vector<std::string>::const_iterator stname = curr_subtracknames.begin();
        stname != curr_subtracknames.end(); ++stname) {
      std::cout << "Processing  " << *stname << std::endl;

      Track<int>::Ptr track(new Track<int>);

      in_file.ReadSubTrack<int>(trackname, *stname, *track);            
      
      NucKDE kde;
      kde.set_out_track_name(trackname);
      kde.set_out_subtrack_name(*stname);
      kde.set_input(track);
      Track<float>::Ptr out_track = kde.output();
      out_file.WriteSubTrack<float>(*out_track);
    }
       
    return 0;
}

