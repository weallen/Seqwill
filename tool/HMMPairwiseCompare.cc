#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

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

void NormalizeTrack(Track<float>* track) {
    float mean = 0.0;
    for (size_t i = 0; i < track->size(); ++i) {
        mean += track->get(i);
    }
    mean /= (float)track->size();
    float temp;
    for (size_t i = 0; i < track->size(); ++i) {
        temp = track->get(i) / mean;
        track->set(i, temp);
    }
}

int main(int argc, char** argv) {
    
    commandArg<std::string> oDataPath("-o", "Out trackfile name");
    commandArg<std::string> iTrack("-i", "In trackfile name");
    commandArg<std::string> tTrackname("-tout", "Out track name");
    commandArg<std::string> t1("-t1", "Track 1");
    commandArg<std::string> t2("-t2", "Track 2");
    commandArg<int> nStates("-n", "Number of states for hmm");
    commandArg<int> wWindows("-w", "Number of windows to merge in output",1);
    
    commandLineParser P(argc, argv);
    P.SetDescription("Compare multiple tracks with an HMM.");
    P.registerArg(oDataPath);
    P.registerArg(iTrack);
    P.registerArg(t1);
    P.registerArg(t2);
    P.registerArg(tTrackname);
    P.registerArg(nStates);
    P.registerArg(wWindows);
    P.parse();
    
    std::string in_trackfile = P.GetStringValueFor(iTrack);
    std::string out_trackfile = P.GetStringValueFor(oDataPath);
    std::string trackname = P.GetStringValueFor(tTrackname);
    std::string track1 = P.GetStringValueFor(t1);
    std::string track2 = P.GetStringValueFor(t2);
    int nstates = P.GetIntValueFor(nStates);
    int win = P.GetIntValueFor(wWindows);
    
    Rng* r = new Rng;
    int num_chrs = 0;
    
    std::cout << "Running HMM with " << nstates << std::endl;
    if (!FileExists(in_trackfile)) {
        std::cerr << "Couldn't find file " << in_trackfile << std::endl;
        delete r;
        return -1;
    }
    
    TrackFile tin(in_trackfile);
    
    if (!tin.HasTrack(track1)) {
        std::cerr << "Couldn't find track " << track1 << std::endl;
        delete r;
        return -1;
    } else {
        std::cout << "Using track " << track1 << std::endl;
    }
    
    if (!tin.HasTrack(track2)) {
        std::cerr << "Couldn't find track " << track2 << std::endl;
        delete r;
        return -1;
    } else {
        std::cout << "Using track " << track2 << std::endl;
    }
    
    std::vector<std::string> chrs = tin.GetSubTrackNames(track1);
    int total_chrlen = 0;
    std::map<std::string, int> chrlens;
    for (size_t i = 0; i < chrs.size(); ++i) {
        TrackMetadata metadata = tin.GetSubTrackMetadata(track1, chrs[i]);
        int len = metadata.GetIntMetadata(chrs[i], 0);
        std::string chrname = chrs[i];
        chrlens[chrname] = len;
        total_chrlen += len;
        std::cout << chrs[i] << " has length " << len << std::endl;
    }
    Track<float>::Ptr alldata(new Track<float>);
    alldata->set_extends(0, total_chrlen);
    
    TrackFile tout(out_trackfile);
    
    GaussHMM cmp(nstates);
    int currpos = 0;
    for (size_t i = 0; i < chrs.size(); ++i) {
        std::string chrname = chrs[i];
        Track<float>* track1_ptr = new Track<float>;
        Track<float>* track2_ptr = new Track<float>;
        
        tin.ReadSubTrack<float>(track1, chrname, *track1_ptr);
        tin.ReadSubTrack<float>(track2, chrname, *track2_ptr);
        
        NormalizeTrack(track1_ptr);
        NormalizeTrack(track2_ptr);
        
        Track<float>* track = new Track<float>;
        track->set_extends(0, floor(track1_ptr->size() / (float)win));
        track->set_resolution(track1_ptr->resolution() * win);
        int len = (int)track1_ptr->size();
        float temp;
        for (size_t i = 0; i < track->size(); ++i) {
            temp = 0.0;
            int start = i*win;
            int stop = std::min((int)(i+1)*win, len);
            for (int j = start; j < stop; ++j) {
                temp += (track1_ptr->get(j) - track2_ptr->get(j));
            }
            track->set(i, temp);
        }
        
        delete track1_ptr;
        delete track2_ptr;
        int currlen = chrlens[chrname];
        for (int j = 0; j < currlen; ++j) {
            alldata->set(currpos + j, track->get(j));
        }
        currpos += currlen + 1;
        delete track;
    }
    
    cmp.set_input(alldata);
    std::vector<GaussDist> dists(cmp.num_states());
    
    Kmeans kmeans(cmp.num_states());
    kmeans.set_track(alldata);
    kmeans.Fit();
    std::vector<double> means = kmeans.means();
    std::sort(means.begin(), means.end());
    for (size_t j = 0; j < cmp.num_states(); ++j) {
        double mean = means[j];
        GaussDist g = GaussDist(mean, 1.0);
        dists[j] = g;
    }
    cmp.set_emit(dists);
    cmp.Init();
    cmp.Compute();
    std::vector<GaussDist> emit = cmp.emit();
    for (size_t i = 0; i < emit.size(); ++i) {
        GaussDist g = emit[i];
        std::cout << "State " << i << " has mean " << g.mean() << std::endl;
    }
    
    Track<int>::Ptr tr = cmp.output();
    
    currpos = 0;
    for (size_t i = 0; i < chrs.size(); ++i) {
        std::cout << "Saving track " << chrs[i] << std::endl;
        Track<int>* currchr = new Track<int>;
        int currlen = chrlens[chrs[i]];
        currchr->set_extends(0, currlen);
        currchr->set_trackname(trackname);
        currchr->set_subtrackname(chrs[i]);
        for (size_t j = 0; j < currlen; ++j) {
            currchr->set(j, tr->get(currpos + j));
        }
        tout.WriteSubTrack<int>(*currchr);
        delete currchr;
        currpos += currlen + 1;
    }
    delete r;
    return 1;
}




