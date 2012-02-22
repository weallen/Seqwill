#include <string>
#include <vector>
#include <iostream>

#include <tbb/task.h>
#include <tbb/concurrent_queue.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/mutex.h>
#include <tbb/atomic.h>



#include "base/CommandLineParser.h"
#include "base/FileUtil.h"
#include "base/Types.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"
#include "analysis/HMM.h"
#include "analysis/Dist.h"
#include "analysis/Random.h"
#include "analysis/Kmeans.h"


int main(int argc, char** argv) {
    
    commandArg<std::string> oDataPath("-o", "Out trackfile name");
    commandArg<std::string> iTrack("-i", "In trackfile name");
    commandArg<std::string> tTrackname("-tout", "Out track name");
    commandArg<std::string> iTrackNames("-tin", "In track names");    
    commandArg<int> nStates("-n", "Number of states for hmm");
    commandArg<int> pThreads("-p", "Number of threads", 1);
    
    commandLineParser P(argc, argv);
    P.SetDescription("Compare multiple tracks with an HMM.");
    P.registerArg(oDataPath);
    P.registerArg(iTrack);
    P.registerCompoundArg(iTrackNames);
    P.registerArg(tTrackname);
    P.registerArg(nStates);
    P.registerArg(pThreads);
    P.parse();
    
    std::string in_trackfile = P.GetStringValueFor(iTrack);
    std::string out_trackfile = P.GetStringValueFor(oDataPath);
    std::string trackname = P.GetStringValueFor(tTrackname);
    std::vector<std::string> in_tracks = P.GetCompoundStringValuesFor(iTrackNames);
    int nstates = P.GetIntValueFor(nStates);
    int threads = P.GetIntValueFor(pThreads);
    
    tbb::task_scheduler_init sched(threads);
    
    tbb::mutex lck;
    Rng* r = new Rng;
    int num_chrs = 0;

    std::cout << "Running HMM with " << nstates << std::endl;
    if (!FileExists(in_trackfile)) {
        std::cerr << "Couldn't find file " << in_trackfile << std::endl;
        delete r;
        return -1;
    }    
    
    TrackFile tin(in_trackfile); 

    
    for (std::vector<std::string>::const_iterator it = in_tracks.begin(); it != in_tracks.end();
         ++it) {
        if (!tin.HasTrack(*it)) {     
            std::cerr << "Couldn't find track " << *it << std::endl;
            delete r;
            return -1;
        } else {
            std::cout << "Using track " << *it << std::endl;
        }
    }
    
    std::vector<std::string> chrs = tin.GetSubTrackNames(in_tracks[0]);
   
    
    TrackFile tout(out_trackfile);    
    for (size_t i = 0; i < chrs.size(); ++i) {
        std::string chrname = chrs[i];
				GaussMultiTrackHMM cmp(nstates);
        cmp.set_out_track_name(trackname);
        cmp.set_out_subtrack_name(chrname);
    
        std::vector<Track<float>::Ptr> chr_tracks;
        for (std::vector<std::string>::iterator it = in_tracks.begin();
             it != in_tracks.end(); ++it) {
						Track<float>::Ptr track(new Track<float>);           
						tin.ReadSubTrack<float>(*it, chrname, *track);
						cmp.add_track(track);
						chr_tracks.push_back(track);
        }
        
        std::vector<std::vector<GaussDist> > dists(cmp.num_states());

        for (int i = 0; i < cmp.num_states(); ++i) {
            std::vector<GaussDist> temp(chr_tracks.size());
            dists[i] = temp;
        }
        
        for (int i = 0; i < chr_tracks.size(); ++i) {
	  Kmeans kmeans(cmp.num_states());
	  kmeans.set_track(chr_tracks[i]);
	  kmeans.Fit() ;
	  std::vector<double> means = kmeans.means();
	  for (size_t j = 0; j < cmp.num_states(); ++j) {
	    double mean = means[j];
	    GaussDist g = GaussDist(mean, 1.0);
	    std::vector<GaussDist> temp = dists[j];
	    temp[i] =  g;
	    dists[j] = temp;
	  }
        }
        cmp.set_emit(dists);
        cmp.Init();
        cmp.Compute();
	Track<int>::Ptr tr = cmp.output();
	std::cout << "Saving track " << tr->subtrackname() << std::endl;
	tout.WriteSubTrack<int>(*tr);
    }

     
    delete r;
    return 1;
}




