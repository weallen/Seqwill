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
#include "analysis/Chromosome.h"

using namespace std;
namespace {
class HDFTrackTest : public ::testing::Test {
  protected:
    HDFTrackTest() {
    }

    virtual ~HDFTrackTest() {
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
  

};
TEST_F(HDFTrackTest, WriteTrackTest) {
}


}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
