#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include "analysis/Genome.h"
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/StringUtil.h"
#include "base/WIG.h"

#include "common/Track.h"


using namespace std;

int main(int argc, char** argv) {
    commandArg<string> iWig("-i", "WIG file");
    commandArg<string> oDataPath("-o", "trackfile name"); 
    commandArg<string> tTrack("-t", "track name");
    commandArg<string> gGenome("-g", "genome info");
    commandArg<int> nSize("-n", "window size");
    commandLineParser P(argc, argv);
    P.SetDescription("Load wig files into HDF5.");
    P.registerArg(iWig);
    P.registerArg(oDataPath);
    P.registerArg(gGenome);
    P.registerArg(nSize);
    P.registerArg(tTrack);
    P.parse();
    
    string wig = P.GetStringValueFor(iWig);
    string dataPath = P.GetStringValueFor(oDataPath);    
    string genomePath = P.GetStringValueFor(gGenome);
    int windowSize = P.GetIntValueFor(nSize);
    string track = P.GetStringValueFor(tTrack);

    GenomeData d;
    d.set_track_name(track);
    GenomeInfo g;
    LoadGenomeInfoFromChr(genomePath, std::string("mm9"), &g);
    d.Init(dataPath, g);
    d.SaveTrackFromWIG(wig, windowSize);
    return 0;
}

