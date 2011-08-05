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
#include "io/BamIO.h"

using namespace std;
namespace {
class BamIOTest : public ::testing::Test {
  protected:
    BamIOTest() {
      b_.Init(std::string("/Users/admin/Documents/omp_nuc_060111.bam"));
    }

    virtual ~BamIOTest() {
    }

    virtual void SetUp() {

    }
    virtual void TearDown() {
    }  
  
  BamIO b_;
};



TEST_F(BamIOTest, SingleReadTest) {
  SingleReadFactory* reads = b_.LoadChrSingleReads(std::string("chr6"));
  int count = 0;
  for (SingleReadFactory::iterator i = reads->begin(); i != reads->end(); ++i) {
    count++;
  }
  ASSERT_EQ(count, 1337874);
  delete reads;
  count = 0;
  reads = b_.LoadAllSingleReads();
  for (SingleReadFactory::iterator i = reads->begin(); i != reads->end(); ++i) {
    count++;
  }

  std::cerr << count << std::endl;
}


}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
