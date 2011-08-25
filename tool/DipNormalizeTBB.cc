#include <string>
#include <vector>
#include <iostream>

#include <tbb/task.h>
#include <tbb/concurrent_queue.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/mutex.h>
#include <tbb/atomic.h>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

#include "base/CommandLineParser.h"
#include "base/FileUtil.h"
#include "base/Types.h"
#include "io/BamIO.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"
#include "analysis/MedipNormalize.h"


tbb::concurrent_queue<Track<float>::Ptr> result_queue;


class NormalizeTask : public tbb::task
{
public:
  NormalizeTask(int extend, int res, const std::string& trackname, 
		const std::string& chrname, TrackFile* chrio,
		BamIO* b, tbb::mutex& lck)
    : chrname_(chrname)
    , b_(b)
    , lck_(lck)
    , chrio_(chrio)
  {
    if (extend > -1) {
      norm_.set_frag_len(extend);
    }
    norm_.set_resolution(res);
    norm_.set_out_track_name(trackname);
    norm_.set_out_subtrack_name(chrname);
  }

  tbb::task* execute()
  {
    Track<unsigned char>::Ptr chr(new Track<unsigned char>);   

    lck_.lock();
    chrio_->ReadSubTrack<unsigned char>(std::string("mm9"), chrname_, *chr);
    SingleReadFactory* reads = b_->LoadChrSingleReads(chr->subtrackname());                 
    lck_.unlock();
    norm_.set_reads(reads);
    norm_.set_chr(chr);
    norm_.Compute();
    result_queue.push(norm_.output());
    delete reads;
    return NULL;
  }

private:
  std::string chrname_;
  BamIO* b_;
  tbb::mutex& lck_;
  TrackFile* chrio_;
  MedipNormalize norm_;
};

int main(int argc, char** argv) {
  
  commandArg<std::string> gGenome("-g", "Genome information");
  commandArg<std::string> bBam("-b", "In BAM file");
  commandArg<std::string> oDataPath("-o", "Trackfile name");
  commandArg<std::string> tTrackname("-t", "Out track name");
  commandArg<int> nRes("-n", "Output resolution", 100); 
  commandArg<int> eExtend("-e", "Extend length if not paired end", -1);
  commandArg<int> pThreads("-p", "Number of threads");
  
  commandLineParser P(argc, argv);
  P.SetDescription("Normalize BAM of me or hmeDIP reads by CpG content.");
  P.registerArg(gGenome);
  P.registerArg(bBam);
  P.registerArg(oDataPath);
  P.registerArg(nRes);
  P.registerArg(eExtend);
  P.registerArg(tTrackname);
  P.registerArg(pThreads);
  P.parse();

  std::string genome = P.GetStringValueFor(gGenome);
  std::string bamname = P.GetStringValueFor(bBam);
  std::string trackfilename = P.GetStringValueFor(oDataPath);
  std::string trackname = P.GetStringValueFor(tTrackname);
  int res = P.GetIntValueFor(nRes);
  int extend = P.GetIntValueFor(eExtend);
  int threads = P.GetIntValueFor(pThreads);

  tbb::task_scheduler_init sched(threads);

  GenomeInfo g;
  BamIO* b = new BamIO(bamname); 

  tbb::mutex lck;

  if (!FileExists(bamname)) {
    std::cerr << "Couldn't find bam file " << bamname << std::endl;
    return -1;
  }
   
  if (!FileExists(genome)) {
    std::cerr << "Couldn't find genome file " << genome << std::endl;
    return -1;
  }

  LoadGenomeInfoFromChr(genome, std::string("mm9"), &g);

  BamTools::RefVector refs = b->reader()->GetReferenceData();  
  std::vector<std::string> chrnames = g.chr_names();

  std::vector<std::string> shared_chrs;
  for (BamTools::RefVector::iterator it = refs.begin(); it != refs.end();
       ++it) {
    if (std::find(chrnames.begin(), chrnames.end(), it->RefName) != chrnames.end()) {     
      if (it->RefName.find("random") == std::string::npos
          && it->RefName.find("NT") == std::string::npos
          && it->RefName.find("chrM") == std::string::npos) {
        shared_chrs.push_back(it->RefName);
      }
    }
  }

  TrackFile* chrio = new TrackFile(genome);
  int num_chrs = (int) shared_chrs.size();


  tbb::empty_task& a = *new(tbb::task::allocate_root()) tbb::empty_task();
  a.set_ref_count(1);

  for (size_t i = 0; i < shared_chrs.size(); ++i) {
    std::string chrname = shared_chrs[i];

    NormalizeTask& task = *new(a.allocate_additional_child_of(a)) NormalizeTask(extend, res, trackname, chrname, chrio, b, lck);
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
  delete b;
  delete chrio;
  return 1;
}




