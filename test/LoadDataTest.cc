#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <vector>
#include <string>

#include "gtest/gtest.h"

#include "base/WIG.h"

using namespace std;
namespace {
class LoadDataTest : public ::testing::Test {
  protected:
    LoadDataTest() {
      filename_ = "/Users/admin/src/Seqwill/testdata/moe_d3a_hmc_50_raw.wig";
    }

    virtual ~LoadDataTest() {
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
  
  std::string filename_;
};
      
TEST_F(LoadDataTest, LoadFromWIGTest) {
  std::vector<WIGLine> lines;
  ParseWig(filename_, &lines);
  ASSERT_EQ(lines[0].val, 0.0);
  ASSERT_EQ(lines[100000].val, 0.778793812);
}

TEST_F(LoadDataTest, WIGToTrackTest) {
  
}

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
