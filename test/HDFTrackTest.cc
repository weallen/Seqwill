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
#include "hdf5/HDFTrackReader.h"
#include "hdf5/HDFTrackWriter.h"

#include "base/SVector.h"

using namespace std;
namespace {
class HDFTrackTest : public ::testing::Test {
  protected:
    HDFTrackTest() {
    }

    virtual ~HDFTrackTest() {
      system("/bin/rm -f /tmp/test.h5");
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
  
    HDFTrackReader<float> reader_;
    HDFTrackWriter<float> writer_;

};
TEST_F(HDFTrackTest, WriterCreateTest) {
  ASSERT_EQ(writer_.Open("/tmp/test.h5"), -1);
  writer_.Close();
  ASSERT_EQ(reader_.Open("/tmp/test.h5"), -1);
  reader_.Close();
}

TEST_F(HDFTrackTest, WriteWriteTrackTest) {

}

TEST_F(HDFTrackTest, WriterGetTrackNamesTest) {

}

TEST_F(HDFTrackTest, ReaderGetTrackNamesTest) {
  reader_.Open("/tmp/test.h5");
  std::vector<std::string> subtracknames = reader_.GetSubtrackNames();
//  ASSERT_EQ(
}

TEST_F(HDFTrackTest, ReaderGetTrackLensTest) {
  reader_.Open("test.h5");
  std::vector<int> subtracklens = reader_.GetSubtrackLengths();
 // ASSERT_EQ(subtracklens.size(), 3);
 // ASSERT_EQ(subtracklens[0], 100);
}

} // namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
