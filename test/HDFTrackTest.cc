#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <vector>
#include <string>

#include "gtest/gtest.h"

#include "data/TrackData.h"

#include "base/SVector.h"
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
TEST_F(HDFTrackTest, ChromosomeTest) {
  SaveChrFromFASTA("out.h5", "/media/Storage/user/data/genomedata/fasta/chr1.fa", "mm9");
}

TEST_F(HDFTrackTest, WriterCreateTest) {
}

TEST_F(HDFTrackTest, WriteWriteTrackTest) {

}

TEST_F(HDFTrackTest, WriterGetTrackNamesTest) {

}

TEST_F(HDFTrackTest, ReaderGetTrackNamesTest) {
//  ASSERT_EQ(
}

TEST_F(HDFTrackTest, ReaderGetTrackLensTest) {
 // ASSERT_EQ(subtracklens.size(), 3);
 // ASSERT_EQ(subtracklens[0], 100);
}

} // namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
