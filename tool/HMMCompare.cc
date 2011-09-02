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


tbb::concurrent_queue<Track<int>::Ptr> result_queue;

class CompareTask : public tbb::task
{
public:
    CompareTask(const std::vector<std::string>& inputs, 
                const std::string& trackname,
                const std::string& chrname, 
                TrackFile* tio, tbb::mutex& lck)
    : trackname_(trackname)
    , chrname_(chrname)
    , lck_(lck)
    {
        cmp_.set_out_track_name(trackname);
        cmp_.set_out_subtrack_name(chrname);
    }
    
    tbb::task* execute()
    {
        Track<float>::Ptr track(new Track<float>);   
        
        lck_.lock();
        tio_->ReadSubTrack<float>(trackname_, chrname_, *track);
        lck_.unlock();
        cmp_.Compute();
        result_queue.push(cmp_.output());
        return NULL;
    }
    
private:
    std::string trackname_;
    std::string chrname_;
    tbb::mutex& lck_;
    TrackFile* tio_;
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
    P.registerArg(pThreads);
    P.parse();
    
    std::string in_trackfile = P.GetStringValueFor(iTrack);
    std::string out_trackfile = P.GetStringValueFor(oDataPath);
    std::string trackname = P.GetStringValueFor(tTrackname);
    std::vector<std::string> in_tracks = P.GetCompoundStringValuesFor(iTrackNames);
    int threads = P.GetIntValueFor(pThreads);
    
    tbb::task_scheduler_init sched(threads);
    
    tbb::mutex lck;
    
    if (!FileExists(in_trackfile)) {
        std::cerr << "Couldn't find file " << in_trackfile << std::endl;
        return -1;
    }    
    
    TrackFile* tin = new TrackFile(in_trackfile); 
    std::vector<std::string> tnames = tin->GetTrackNames();
    
    for (std::vector<std::string>::const_iterator it = in_tracks.begin(); it != in_tracks.end();
         ++it) {
        if (!tin->HasTrack(*it)) {     
            std::cerr << "Couldn't fine track " << *it << std::endl;
            return -1;
        }
    }
    
    std::vector<std::string> chrs = tin->GetSubTrackNames(tnames[0]);
    
    tbb::empty_task& a = *new(tbb::task::allocate_root()) tbb::empty_task();
    a.set_ref_count(1);
    
    for (size_t i = 0; i < chrs.size(); ++i) {
        std::string chrname = chrs[i];
        
        CompareTask& task = *new(a.allocate_additional_child_of(a)) CompareTask(in_tracks, trackname, chrname, tio, lck);
        a.spawn(task);
    }
    std::cout << "All tasks enqueued." << std::endl;
    
    TrackFile tio(trackfilename);
    
    int chrs_done = 0;
    while(1) {
        if (!result_queue.empty()) {
            Track<float>::Ptr tr;
            result_queue.try_pop(tr);
            std::cout << "Saving track " << tr->subtrackname() << std::endl;
            tio.WriteSubTrack<float>(*tr);
            chrs_done++;
        }
        if (chrs_done == num_chrs) {
            break;
        }
    }
    
    a.wait_for_all();
    a.destroy(a);
    
    std::cout << "Cleaning up" << std::endl;
    
    delete tin;
    return 1;
}




