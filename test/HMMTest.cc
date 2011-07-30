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
      track2_ = Track<float>::Ptr(new Track<float>);
      h_ = new GaussHMM(2);
      h2_ = new BernHMM(3);
      TrackFile f;
      TrackFile f2;
      f.Open(std::string("/Users/admin/Documents/test_hmm.trk"));
      f.ReadSubTrack<float>(std::string("testdata"), std::string("test1"), *track_);

      f2.Open(std::string("/Users/admin/Documents/test.trk"));
      f2.ReadSubTrack<float>(std::string("moe_d3a_hmc_raw"), std::string("chr1"), *track2_);          
      h_->set_input(track_);
      
      HMM::MatrixType trans_prior = HMM::MatrixType::Constant(3, 3, 1/3);
      HMM::VectorType init = HMM::VectorType::Constant(3, 1/3);
      //h2_->set_transition(trans_prior);
      //h2_->set_trans_prior(trans_prior);
      //h2_->set_init_probs(init);
      //h2_->set_init_probs_prior(init);
      
      std::vector<GaussDist> g;
      g.push_back(GaussDist(5.0, 1.0));
      g.push_back(GaussDist(5.0, 1.0));
      h_->set_emit(g);
      h_->Init();
            
      int num = 0;
      for (int i = 0; i < track2_->size(); ++i) {
        if (track2_->get(i) > 0.3) {
          track2_->set(i, 1.0);
          num++;
        } else {
          track2_->set(i,0.0);
        }
      }
      std::cerr << "NUM > 0.3: "<< num << " out of " << track2_->size() << std::endl;
      h2_->set_input(track2_);

      std::vector<BernDist> g2;
      HMM::VectorType means = HMM::VectorType::Random(3);

      for (int i = 0; i < 3; ++i) {
        g2.push_back(BernDist(means(i)));
      }
      h2_->set_emit(g2);
      h2_->Init();
    }
    
    virtual ~HMMTest() {
      delete h_;
      delete h2_;
    }
    
    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
    
    Track<float>::Ptr track_;
    Track<float>::Ptr track2_;
    GaussHMM* h_;
    BernHMM* h2_;
  };
  
  
  TEST_F(HMMTest, FitRealDataTest) {
    std::cerr << h2_->transition() << std::endl;
    h2_->FitEM();
    HMM::StateVectorType path;
    h2_->Decode(path);
    std::cerr << h2_->transition() << std::endl;
    int num = 0;
    for (int i = 0; i < path.cols(); ++i) {
      if (path(i) > 0) {
        num++;
      }
    }
    std::cerr << "NUM GREATER THAN 0 " << num << std::endl;
//    std::cerr << (path == 0).count() << std::endl;
//    std::cerr << (path == 1).count() << std::endl;
//    std::cerr << (path == 2).count() << std::endl;
  }

  TEST_F(HMMTest, FitEMTest) {
    h_->FitEM();
    std::vector<GaussDist> e = h_->emit();
    std::cerr << e[0].mean() << " , " << e[0].stddev() << std::endl;
    std::cerr << e[1].mean() << " , " << e[1].stddev() << std::endl;
    HMM::StateVectorType path;
    h_->Decode(path);
    ASSERT_EQ((path == 0).count(),12520 );
    ASSERT_EQ((path == 1).count(), 87480);
    std::cerr << h_->transition() << std::endl;
  }
  
}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
