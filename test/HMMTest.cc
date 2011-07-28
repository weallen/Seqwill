#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <vector>
#include <string>

#include "gtest/gtest.h"

#include "common/Track.h"
#include "io/TrackIO.h"

#include "analysis/Dist.h"
#include "analysis/HMM.h"

using namespace std;
namespace {
class HMMTest : public ::testing::Test {
  protected:
    HMMTest() 
  {
    track_ = Track<float>::Ptr(new Track<float>);
      h_ = new GaussHMM();
      TrackFile f;
      f.Open(std::string("/Users/admin/Documents/test_hmm.trk"));
      f.ReadSubTrack(std::string("testdata"), std::string("test1"), *track_);
      h_->set_input(track_);
      h_->set_num_states(2);
      
      HMM::MatrixType trans = HMM::MatrixType::Constant(2, 2, 0.5);
      HMM::VectorType init = HMM::VectorType::Constant(2, 0.5);

      h_->set_transition(trans);
      h_->set_init_probs(init);
    }

    virtual ~HMMTest() {
      delete h_;
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
  
  Track<float>::Ptr track_;
  GaussHMM* h_;
};
  

TEST_F(HMMTest, LogProbTest) {
  HMM::MatrixType softev = HMM::MatrixType::Random(2, 10);
  std::cerr << h_->LogProb(softev) << std::endl;
}

TEST_F(HMMTest, FwdBackTest) {
  h_->FitEM();
}

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
