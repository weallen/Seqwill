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
    HMMTest() {
      TrackFile f;
      f.Open(std::string("/Users/admin/Documents/test_hmm.trk"));
      f.ReadSubTrack(std::string("testdata"), std::string("test1"), *track_);
      h_.set_num_states(2);
    }

    virtual ~HMMTest() {
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
  
  Track<float>::Ptr track_;
  HMM h_;
};
  
TEST_F(HMMTest, InitHMMTest) {
  
}


}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
