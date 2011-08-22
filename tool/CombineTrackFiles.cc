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
#include "base/WIG.h"

#include "common/Track.h"
#include "io/TrackIO.h"

using namespace std;

int main(int argc, char** argv) {
    commandArg<string> oDataPath("-o", "Out trackfile name"); 
    commandArg<string> tTrackfiles("-t", "In trackfiles");

    commandLineParser P(argc, argv);
    P.SetDescription("Combine tracks into one HDF5.");
    P.registerArg(oDataPath);
    P.registerCompoundArg(tTrackfiles);
    P.parse();
    
    string dataPath = P.GetStringValueFor(oDataPath);    
    vector<string> trackfiles = P.GetCompoundStringValuesFor(tTrackfiles);

    TrackFile tio(dataPath);
    for (size_t i = 0; i < trackfiles.size(); ++i) {
        std::string curr_track = trackfiles[i]; 
        TrackFile curr_file(curr_track);
        std::vector<std::string> curr_tracknames = curr_file.GetTrackNames();
        for (std::vector<std::string>::const_iterator tname = curr_tracknames.begin();
            tname != curr_tracknames.end(); ++tname) {
            std::vector<std::string> curr_subtracknames = curr_file.GetSubTrackNames(*tname);
            for (std::vector<std::string>::const_iterator stname = curr_subtracknames.begin();
                stname != curr_subtracknames.end(); ++stname) { 
                Track<float>::Ptr track(new Track<float>);
                curr_file.ReadSubTrack<float>(*tname, *stname, *track);
                tio.WriteSubTrack<float>(*track);
            }
        }
    }

    return 0;
}

