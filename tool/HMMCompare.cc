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

tbb::concurrent_queue<Track<int>::Ptr> result_queue;

class CompareTask : public tbb::task
{
public:
    CompareTask(std::vector<std::string>& inputs, 
                std::string& trackname,
                std::string& chrname, 
                TrackFile* tio, tbb::mutex& lck,
                Rng* r, int num_states)
    : inputs_(inputs)
    , trackname_(trackname)
    , chrname_(chrname)
    , lck_(lck)
    , tio_(tio)
    , r_(r)
    {
        cmp_.set_out_track_name(trackname);
        cmp_.set_out_subtrack_name(chrname);
        cmp_.set_num_states(num_states);
    }
    
    tbb::task* execute()
    {
        for (std::vector<std::string>::iterator it = inputs_.begin();
             it != inputs_.end(); ++it) {
            Track<float>::Ptr track(new Track<float>);           
            lck_.lock();
            tio_->ReadSubTrack<float>(*it, chrname_, *track);
            lck_.unlock();
            cmp_.add_track(track);
        }
        
        std::vector<std::vector<GaussDist> > dists(cmp_.num_states());

        for (size_t i = 0; i < cmp_.num_states(); ++i) {
            for (int j = 0; j < inputs_.size(); ++j) {
                lck_.lock();
               double mean = 10*abs(gsl_rng_uniform(r_->rng()));
                lck_.unlock();
               std::cout << mean << std::endl;
                dists[i].push_back(GaussDist(mean,1.0));
            }
        }
        cmp_.set_emit(dists);
        cmp_.Init();
        cmp_.Compute();
        result_queue.push(cmp_.output());
        return NULL;
    }
    
private:
    std::vector<std::string>& inputs_;
    std::string& trackname_;
    std::string& chrname_;
    tbb::mutex& lck_;
    TrackFile* tio_;
    Rng* r_;
    GaussMultiTrackHMM cmp_;
};

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
    
    tbb::empty_task& a = *new(tbb::task::allocate_root()) tbb::empty_task();
    a.set_ref_count(1);
    
    
    for (size_t i = 0; i < chrs.size(); ++i) {
        std::string chrname = chrs[i];
        num_chrs++;
        CompareTask& task = *new(a.allocate_additional_child_of(a)) CompareTask(in_tracks, trackname, chrname, &tin, lck, r, nstates);
        a.spawn(task);
    }
    std::cout << "All tasks enqueued." << std::endl;
    
    TrackFile tout(out_trackfile);
    int chrs_done = 0;
    while(1) {
        if (!result_queue.empty()) {
            Track<int>::Ptr tr;
            result_queue.try_pop(tr);
            std::cout << "Saving track " << tr->subtrackname() << std::endl;
            tout.WriteSubTrack<int>(*tr);
            chrs_done++;
        }
        if (chrs_done == num_chrs) {
            break;
        }
    }
    
    a.wait_for_all();
    a.destroy(a);
    
    std::cout << "Cleaning up" << std::endl;
    
    delete r;
    return 1;
}




