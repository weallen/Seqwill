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
      bio_ = new BamIO(std::string("/Users/wea/Desktop/icam_nuc_012.bam"));
        
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
  
/*TEST_F(NucMapTest, PileupTest) {
      NucPileup pileup;
      pileup.set_out_track_name(std::string("test"));
      pileup.set_out_subtrack_name(std::string("chr19"));
      SingleReadFactory* reads = bio_->LoadChrSingleReads(std::string("chr19"));
      pileup.set_reads(reads);
      pileup.set_start(0);
      pileup.set_stop(149517037);
      pileup.Compute();
			std::cerr << "Mean " << pileup.fraglen_mean() << std::endl;
			std::cerr << "Stddev " << pileup.fraglen_var() << std::endl;
      std::cerr << "Got here" << std::endl;
      TrackFile tio("/Users/wea/Desktop/nuc_map_test.trk");
      Track<int>::Ptr track = pileup.output();
      tio.WriteSubTrack<int>(*track);
      delete reads; 
			std::fstream out;
		  out.open("/Users/wea/Desktop/nuc.wig", std::fstream::out);
			out << "fixedStep chrom=chr6 start=0 step=25 span=25" << std::endl;
			for (int i = 0; i < (int)track->size(); i += 25) {
				int sum = 0;
				for (int j = i; j < (i+25); ++j) {
					sum += (float)track->get(j);
				}
				out << sum << std::endl;
			}
  }
  
  
*/
  TEST_F(NucMapTest, KDETest) {
    NucKDE kde;
    Track<int>::Ptr track(new Track<int>);

    TrackFile tio("/Users/wea/Desktop/nuc_map_test.trk");
    tio.ReadSubTrack<int>(std::string("test"), std::string("chr19"), *track);
    int s = 0;
    for (size_t i = 0; i < track->size(); ++i) {
      s += track->get(i);
    }
    std::cout << s << std::endl;
      kde.set_out_track_name("pos_test");
      kde.set_out_subtrack_name("chr19");
    kde.set_input(track);
    kde.Compute();
    Track<float>::Ptr out = kde.output();
	float sum = 0.0;
	for (size_t i = 0; i < out->size(); ++i) {
		sum += out->get(i);
	}
	std::cout << sum << std::endl;
    tio.WriteSubTrack<float>(*out);
  }
  

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
