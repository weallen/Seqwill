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
    commandArg<string> tTrackfiles("-t", "In trackfile");
    commandArg<string> nName("-n", "Track name");
    commandArg<string> wWig("-w", "Out wig file");
    commandLineParser P(argc, argv);
    P.SetDescription("Turn track into WIG.");
    P.registerArg(wWig);
    P.registerArg(nName);
    P.registerArg(tTrackfiles);
    P.parse();
    
    string wig = P.GetStringValueFor(wWig);    
    string trackname = P.GetStringValueFor(nName);
    string trackfile = P.GetStringValueFor(tTrackfiles);

    
    
    TrackFile curr_file(trackfile);
    std::vector<std::string> tnames = curr_file.GetTrackNames();
    if (std::find(tnames.begin(), tnames.end(), trackname) == tnames.end()) {
        std::cerr << "Couldn't find track " << trackname << std::endl;
        return -1;
    }

    std::fstream out;
    out.open((wig + ".wig").c_str(), std::fstream::out);

    std::vector<std::string> curr_subtracknames = curr_file.GetSubTrackNames(trackname);
    for (std::vector<std::string>::const_iterator stname = curr_subtracknames.begin();
        stname != curr_subtracknames.end(); ++stname) { 
        Track<float>::Ptr track(new Track<float>);
        curr_file.ReadSubTrack<float>(trackname, *stname, *track);
        out << "fixedStep chrom=" << *stname  << " start=0 step=" << track->resolution() << std::endl;
        for (size_t i = 0; i < track->stop(); ++i) {
            out << track->get(i) << std::endl;
        }
    }
       
    out.close();
    return 0;
}

