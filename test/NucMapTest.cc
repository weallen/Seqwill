#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "gtest/gtest.h"

#include "common/Track.h"
#include "io/BamIO.h"
#include "analysis/Nucleosome.h"
#include "io/TrackIO.h"

using namespace std;
namespace {
  class NucMapTest : public ::testing::Test {
  protected:
    NucMapTest() {
      bio_ = new BamIO(std::string("/media/storage2/data/bam/sorted_cells/nuc/omp_nuc_060111.bam"));
        
    }
    
    virtual ~NucMapTest() {
      delete bio_;
    }
    
    virtual void SetUp() {
      
    }
    virtual void TearDown() {
    }  
    
    BamIO* bio_;
  };
  
  TEST_F(NucMapTest, PileupTest) {
      NucPileup pileup;
      pileup.set_out_track_name(std::string("test"));
      pileup.set_out_subtrack_name(std::string("chr6"));
      SingleReadFactory* reads = bio_->LoadChrSingleReads(std::string("chr6"));
      pileup.set_reads(reads);
      pileup.Compute();
      TrackFile tio("/home/will/Documents/nuc_map_test.trk");
      tio.WriteSubTrack<int>(pileup.output());
      delete reads;
      
  }

  
  
  TEST_F(NucMapTest, AssignCpGToFragTest) {
  }
  

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
