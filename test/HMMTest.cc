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
#include "analysis/Kmeans.h"

using namespace std;
namespace {
  class HMMTest : public ::testing::Test {
  protected:
    HMMTest() 
    {
      track_ = Track<float>::Ptr(new Track<float>);
      track2_ = Track<float>::Ptr(new Track<float>);
      track3_ = Track<float>::Ptr(new Track<float>);
      track4_ = Track<float>::Ptr(new Track<float>);
      track5_ = Track<float>::Ptr(new Track<float>);
      track6_ = Track<float>::Ptr(new Track<float>);
      h_ = new GaussHMM(2);
      h2_ = new GaussHMM(6);
      TrackFile f;
      TrackFile f2;
      f.Open(std::string("/Users/admin/Documents/test2.trk"));
      f.ReadSubTrack<float>(std::string("testdata"), std::string("test1"), *track_);
      f.ReadSubTrack<float>(std::string("testdata"), std::string("test2"), *track4_);
      f2.Open(std::string("/Users/admin/Documents/test.trk"));
      f2.ReadSubTrack<float>(std::string("moe_d3a_hmc_50_raw"), std::string("chr6"), *track2_);          
      f2.ReadSubTrack<float>(std::string("moe_wt_hmc_50_raw"), std::string("chr6"), *track3_);
      f2.ReadSubTrack<float>(std::string("moe_wt_mc_50_raw"), std::string("chr6"), *track5_);
      f2.ReadSubTrack<float>(std::string("moe_d3a_mc_50_raw"), std::string("chr6"), *track6_);

      h_->set_input(track_);
      
      int num = 0;
      
      for (size_t i = 0; i < track2_->size(); ++i) {
        double temp = abs(track2_->get(i) - track3_->get(i));
        track2_->set(i, temp);
      }
      std::cerr << num << " greater than 0.5, out of " << track2_->size() << std::endl;
      
      h2_->set_input(track2_);
    }
    
    virtual ~HMMTest() {
      delete h_;
      delete h2_;
    }
    
    virtual void SetUp() {
      HMM::MatrixType trans_prior = HMM::MatrixType::Constant(3, 3, 1/3);
      HMM::VectorType init = HMM::VectorType::Constant(3, 1/3);
      //h2_->set_transition(trans_prior);
      //h2_->set_trans_prior(trans_prior);
      //h2_->set_init_probs(init);
      //h2_->set_init_probs_prior(init);
      
      std::vector<GaussDist> g;
      g.push_back(GaussDist(4.5, 1.0));
      g.push_back(GaussDist(5.5, 1.0));
      h_->set_emit(g);
      h_->Init();
            
      std::vector<GaussDist> g2;
      HMM::VectorType means = HMM::VectorType::Random(6);
      
      for (int i = 0; i < 6; ++i) {
        g2.push_back(GaussDist(abs(means(i)),1.0));
      }
      h2_->set_emit(g2);
      h2_->Init();

    }
    virtual void TearDown() {
    }
    
    Track<float>::Ptr track_;
    Track<float>::Ptr track2_;
    Track<float>::Ptr track3_;
    Track<float>::Ptr track4_;
    Track<float>::Ptr track5_;
    Track<float>::Ptr track6_;
    GaussHMM* h_;
    GaussHMM* h2_;
  };

  TEST_F(HMMTest, FitRealDataTest) {
    GaussMultiTrackHMM hmm(5);
    hmm.add_track(track2_);
    hmm.add_track(track3_);
    hmm.add_track(track5_);
    hmm.add_track(track6_);
    std::vector<std::vector<GaussDist> > g(5);
    std::vector<GaussDist> temp(4);
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(5, 4);
    for (int i = 0; i < 5; ++i) {
      temp[0] = GaussDist(m(i, 0), 1.0);
      temp[1] = GaussDist(m(i, 1), 1.);
      g[i] = temp;
    }
    hmm.set_emit(g);
    hmm.Init();
      hmm.FitEM();
    //h2_->FitBlockedGibbs();
    HMM::StateVectorType path;
    hmm.Decode(path);
    int num = 0;
    for (int i = 0; i < path.cols(); ++i) {
      if (path(i) > 0) {
        num++;
      }
    }
    std::fstream f("/Users/admin/Documents/output_gauss.wig",std::fstream::out);
    f << "fixedStep  chrom=chr6  start=0  step=50\n";
    for (int i = 0; i < path.size(); ++i) {
      f << path(i) << std::endl;
    }
  }
    TEST_F(HMMTest, FitRealDataTestMV) {
    //std::cerr << h2_->transition() << std::endl;
      Kmeans k(10);

      k.add_track(track2_);
      k.add_track(track3_);

      //std::vector<Eigen::VectorXd> means = k.Fit();
      MVGaussMultiTrackHMM hmm(10);
      hmm.add_track(track2_);
      hmm.add_track(track3_);
      std::vector<MVGaussDist> g;
      Eigen::MatrixXd means = Eigen::MatrixXd::Random(10,4);
      for (int i = 0; i < 10; ++i) {
	std::cerr << "Mean " << i << " = " << means.row(i).cwiseAbs() << std::endl;
	g.push_back(MVGaussDist(means.row(i).cwiseAbs(), Eigen::MatrixXd::Identity(2,2)));
      }
      hmm.set_emit(g);
      hmm.Init();
      hmm.FitEM();
    //h2_->FitBlockedGibbs();
    HMM::StateVectorType path;
    hmm.Decode(path);
    int num = 0;
    for (int i = 0; i < path.cols(); ++i) {
      if (path(i) > 0) {
        num++;
      }
    }
    std::fstream f("/Users/admin/Documents/output.wig",std::fstream::out);
    f << "fixedStep  chrom=chr6  start=0  step=50\n";
    for (int i = 0; i < path.size(); ++i) {
      f << path(i) << std::endl;
    }
    std::cerr << "NUM GREATER THAN 0 " << num << std::endl;
    std::cerr << (path == 0).count() << std::endl;
    std::cerr << (path == 1).count() << std::endl;
    std::cerr << (path == 2).count() << std::endl;
  }


  TEST_F(HMMTest, MultiTrackTest) {
    MVGaussMultiTrackHMM hmm(2);
    hmm.add_track(track_);
    hmm.add_track(track4_);
    hmm.Init();
    std::vector<MVGaussDist> g(2);
    HMM::VectorType means = HMM::VectorType::Random(2);

    for (int i = 0; i < 2; ++i) {
      g.push_back(MVGaussDist(Eigen::VectorXd::Random(2).cwiseAbs(),Eigen::MatrixXd::Identity(2,2)));
      std::cerr << g[i].mean() << std::endl;
    }
    hmm.set_emit(g);
    hmm.FitEM();
    HMM::StateVectorType path;
    h_->Decode(path);
    ASSERT_EQ((path == 0).count(), 87480);
    ASSERT_EQ((path == 1).count(), 12520 );
    std::cerr << h_->transition() << std::endl;
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

  

  TEST_F(HMMTest, FitGibbsTest) {
    h_->FitBlockedGibbs();
    HMM::StateVectorType path;
    h_->Decode(path);
    std::cerr << (path == 0).count() << " in state 0" << std::endl;
    std::cerr << (path == 1).count() << " in state 1" << std::endl;
  }
  

  
}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
