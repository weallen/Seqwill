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

#include "io/TrackIO.h"
#include "common/Track.h"
#include "analysis/Kmeans.h"

using namespace std;
namespace {
  class KmeansTest : public ::testing::Test {
  protected:
    KmeansTest() {
      fname_ = std::string("/media/storage2/data/h5/rpm_avg/all_dips_rpm_avg.trk");
        
    }
    
    virtual ~KmeansTest() {
    }
    
    virtual void SetUp() {
      
    }
    virtual void TearDown() {
    }  
    
    std::string fname_;
  };
  
  TEST_F(KmeansTest, KTest) {
    Kmeans kmeans(3);
    Track<float>::Ptr track(new Track<float>);
    TrackFile tio(fname_);
    tio.ReadSubTrack<float>(std::string("icam_hmc"), std::string("chr19"), *track);
    kmeans.set_track(track);
    kmeans.Fit();
    std::vector<double> means = kmeans.means();
    for (int i = 0; i < 3; ++i) {
     std::cout << means[i] << std::endl;
    }
  }
  

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
