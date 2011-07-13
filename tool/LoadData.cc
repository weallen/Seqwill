#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/StringUtil.h"
#include "base/WIG.h"

#include "common/Track.h"


using namespace std;

int main(int argc, char** argv) {
    commandArg<string> iWig("-i", "WIG file");
    commandArg<string> oDataPath("-o", "trackname"); 
    commandLineParser P(argc, argv);
    P.SetDescription("Load wig files into HDF5.");
    P.registerArg(iWig);
    P.registerArg(oDataPath);
    P.parse();
    
    std::string wig = P.GetStringValueFor(iWig);
    std::string dataPath = P.GetStringValueFor(oDataPath);
        
    return 0;
}

