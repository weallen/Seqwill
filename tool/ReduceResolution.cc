#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include "base/CommandLineParser.h"
#include "base/FileUtil.h"
#include "base/Types.h"
#include "io/TrackIO.h"

int main(int argc, char** argv) {

    commandArg<std::string> iDataPath("-i", "In trackfile");
    commandArg<std::string> tTrack("-t", "In trackname");
    commandArg<std::string> oDataPath("-o", "Trackfile name");
    commandArg<int> nFactor("-n", "Factor to reduce resolution by");
    commandLineParser P(argc, argv);
    P.SetDescription("Reduce resolution by summing together windows");

    P.registerArg(oDataPath);
    P.registerArg(iDataPath);
    P.registerArg(nFactor);
    P.registerArg(tTrack);
    P.parse();

    std::string outtrackfilename = P.GetStringValueFor(oDataPath);
    std::string intrackfilename = P.GetStringValueFor(iDataPath);
    std::string trackname = P.GetStringValueFor(tTrack);
    int res = P.GetIntValueFor(nFactor);

    if (!FileExists(intrackfilename)) {
        std::cerr << "Couldn't find track file " << intrackfilename << std::endl;
        return -1;
    }


    TrackFile tin(intrackfilename);
    TrackFile tout(outtrackfilename);

    std::vector<std::string> tnames = tin.GetSubTrackNames(trackname);
    for (std::vector<std::string>::iterator it = tnames.begin(); it != tnames.end(); ++it) {
        std::string stname = *it;
        Track<float>::Ptr track(new Track<float>);
        tin.ReadSubTrack(trackname, stname, *track);

        Track<float>::Ptr newtrack(new Track<float>);
        newtrack->set_extends(0, static_cast<int>(ceil(track->stop()/static_cast<float>(res))));
        newtrack->set_trackname(track->trackname());
        newtrack->set_subtrackname(stname);
        newtrack->set_resolution(track->resolution() * res);
        for (int i = 0; i < newtrack->stop(); ++i) {
            newtrack->set(i, 0.0);
        }
        for (int i = 0; i < newtrack->stop(); ++i) {
            for (int j = (i * res); j < ((i+1) * res); ++j) {
                if (j < track->stop()) {
                    newtrack->set(i, newtrack->get(i) + track->get(j));
                }
            }
        }
        tout.WriteSubTrack(*newtrack);
    }
    tin.Close();
    tout.Close();
    return 1;
}
