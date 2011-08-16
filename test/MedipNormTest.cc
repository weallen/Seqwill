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

#include "gtest/gtest.h"

#include "common/Track.h"
#include "io/TrackIO.h"
#include "analysis/MedipNormalize.h"

using namespace std;
namespace {
class MedipNormTest : public ::testing::Test {
  protected:
    MedipNormTest() {
      srand(time(NULL));
      chr6_ = Track<unsigned char>::Ptr(new Track<unsigned char>);
      tio_.Open("/Users/admin/Documents/test.trk");
      chr_.Open("/Users/admin/Documents/chromosomes.trk");
      chr_.ReadSubTrack(std::string("mm9"), std::string("chr6"), *chr6_);
    }

    virtual ~MedipNormTest() {
    }

    virtual void SetUp() {

    }
    virtual void TearDown() {
    }  
  TrackFile tio_;
  TrackFile chr_;
  Track<unsigned char>::Ptr chr6_;
};

  TEST_F(MedipNormTest, CouplingFactorTest) {
    MedipNormalize norm;
    Track<float>::Ptr t = Track<float>::Ptr(new Track<float>);
    tio_.ReadSubTrack(std::string("moe_wt_hmc_50_raw"), std::string("chr6"), *t);
    norm.set_input(t);
    norm.set_chr(chr6_);
    norm.set_frag_len(300);
    Eigen::VectorXd v;
    norm.Compute();
    v = norm.coupling();
    std::fstream f;
    f.open("/Users/admin/Documents/chr6_couple.txt", std::fstream::out);
    for (int i = 0; i < v.size(); ++i) {
      f << v(i) << std::endl;
    }
  }

  TEST_F(MedipNormTest, ReadCompoundTrackTest) {
  }

TEST_F(MedipNormTest, WriteTrackTest) {
}

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
