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
#include "analysis/MedipNormalize.h"
#include "io/TrackIO.h"

using namespace std;
namespace {
  class MedipNormTest : public ::testing::Test {
  protected:
    MedipNormTest() {
      bio_ = new BamIO(std::string("/Users/admin/Documents/d3a_wt_mc.bam"));
      chr6_ = Track<unsigned char>::Ptr(new Track<unsigned char>);
      chr_.Open("/Users/admin/Documents/chromosomes.trk");
      chr_.ReadSubTrack(std::string("mm9"), std::string("chr6"), *chr6_);
      norm_.set_chr(chr6_);
      norm_.set_bam(bio_);
    }
    
    virtual ~MedipNormTest() {
      delete bio_;
    }
    
    virtual void SetUp() {
      
    }
    virtual void TearDown() {
    }  
    
    MedipNormalize norm_;
    TrackFile chr_;
    Track<unsigned char>::Ptr chr6_;
    BamIO* bio_;
  };
  
  TEST_F(MedipNormTest, ComputeTest) {
    TrackFile tio(std::string());
    
    norm_.set_out_track_name(std::string("d3a_wt_mc"));
    norm_.set_out_subtrack_name(std::string("chr6"));
    norm_.Compute();
    Track<float>::Ptr t = norm_.output();
    std::fstream outfile;
    outfile.open("/Users/admin/Documents/norm_mc_test.wig", std::fstream::out);
    outfile << "fixedStep chrom=chr6 start=0 step=100" << std::endl;
    for (size_t i = 0; i < t->size(); ++i) {
      outfile << t->get(i) << std::endl;
    }    
    outfile.close();
  }
  
  TEST_F(MedipNormTest, AssignCpGToFragTest) {
    norm_.FindCpG();
    ASSERT_EQ(norm_.num_cpgs(), 1159498);
    norm_.ReadsToFrags();
    ASSERT_EQ(norm_.frags().size(), 930662);    
    norm_.AssignCpGToFrags();
    std::vector<Frag> frags = norm_.frags();
    int total = 0;
    for (size_t i = 0; i < frags.size(); ++i) {
      total += frags[i].num_cpgs;
    }
    ASSERT_EQ(total, 3748724);    
  }
  

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
