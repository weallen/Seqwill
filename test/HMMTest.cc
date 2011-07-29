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
      h_ = new GaussHMM(2);
      TrackFile f;
      f.Open(std::string("/Users/admin/Documents/test_hmm.trk"));
      f.ReadSubTrack(std::string("testdata"), std::string("test1"), *track_);
      h_->set_input(track_);
      
      HMM::MatrixType trans_prior = HMM::MatrixType::Constant(2, 2, 0.25);
      HMM::VectorType init = HMM::VectorType::Constant(2, 0.5);

      h_->set_trans_prior(trans_prior);
      h_->set_init_probs_prior(init);
    std::vector<GaussDist> g;
    g.push_back(GaussDist(10.0, 1.0));
    g.push_back(GaussDist(0.0, 1.0));
    h_->set_emit(g);

      h_->Init();
    
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
  
TEST_F(HMMTest, FitEMTest) {
  h_->FitEM();
  std::vector<GaussDist> e = h_->emit();
  std::cerr << e[0].mean() << " , " << e[0].stddev() << std::endl;
  std::cerr << e[1].mean() << " , " << e[1].stddev() << std::endl;
}

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
