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
#include "base/BEDelement.h"

using namespace std;

int main(int argc, char** argv) {
    commandArg<string> tTrackfile("-i", "In trackfile");
    commandArg<string> oFile("-o", "Out BED file");
    commandArg<string> tTrack("-t", "Track name");

    commandLineParser P(argc, argv);
    P.SetDescription("Find nucleosome positioning stringency.");
    P.registerArg(tTrackfile);
    P.registerArg(oFile);
    P.registerArg(tTrack);
    P.parse();
    

    string trackname = P.GetStringValueFor(tTrack);
    string in_trackfile = P.GetStringValueFor(tTrackfile);
    string out_file = P.GetStringValueFor(oFile);
    
    TrackFile in_file(in_trackfile);

    std::fstream out;
    out.open(out_file.c_str(), std::fstream::out);

    std::vector<std::string> tnames = in_file.GetTrackNames();
    if (std::find(tnames.begin(), tnames.end(), trackname) == tnames.end()) {
        std::cerr << "Couldn't find track " << trackname << std::endl;
        return -1;
    }
    
    
    std::vector<std::string> curr_subtracknames = in_file.GetSubTrackNames(trackname);
    for (std::vector<std::string>::const_iterator stname = curr_subtracknames.begin();
        stname != curr_subtracknames.end(); ++stname) {
      std::cout << "Processing  " << *stname << std::endl;

      Track<float>::Ptr track(new Track<float>);

      in_file.ReadSubTrack<float>(trackname, *stname, *track);            
      
      NucPositioner p;
      p.set_chrname(*stname);
      p.set_input(track);
      p.FindNucs();
      std::vector<BEDelement> elems = p.NucPosAsBED();
      for (std::vector<BEDelement>::iterator it = elems.begin();
	   it != elems.end(); ++it) {
	out << *it << std::endl;
      }     
    }

    out.close();
    return 0;
}

