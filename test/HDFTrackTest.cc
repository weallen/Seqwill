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
#include "io/TrackIO.h"
#include "analysis/Chromosome.h"
#include "io/Traits.h"
#include "base/Types.h"

using namespace std;
namespace {
class HDFTrackTest : public ::testing::Test {
  protected:
    HDFTrackTest() {
      srand(time(NULL));
      tio_.Open("/Users/admin/Documents/test_hdf.trk");
    }

    virtual ~HDFTrackTest() {
    }

    virtual void SetUp() {

    }
    virtual void TearDown() {
    }  
  TrackFile tio_;
};

  TEST_F(HDFTrackTest, WriteCompoundTrackTest) {
    std::string s = DataTypeTraits<PlusMinusDataFloat>::Name();
    ASSERT_EQ(s, "PlusMinusDataFloat");
    Track<PlusMinusDataFloat>::Ptr t(new Track<PlusMinusDataFloat>());
    t->set_extends(0,100);
    t->set_resolution(1);
    t->set_trackname("test");
    t->set_subtrackname("test_compound");
    PlusMinusDataFloat d;
    for (int i = 0; i < 100; ++i) {
      d.plus = 1.0;
      d.minus = 0.0;
      t->set(i, d);
    }
    ASSERT_TRUE(tio_.WriteSubTrack<PlusMinusDataFloat>(*t));
  }

  TEST_F(HDFTrackTest, ReadCompoundTrackTest) {
    Track<PlusMinusDataFloat>::Ptr t(new Track<PlusMinusDataFloat>());
    tio_.ReadSubTrack<PlusMinusDataFloat>(std::string("test"), std::string("test_compound"),
					  *t);
    PlusMinusDataFloat d = t->get(0);
    ASSERT_EQ(d.plus, 1.0);
    ASSERT_EQ(d.minus, 0.0);
  }

TEST_F(HDFTrackTest, WriteTrackTest) {
  Track<float>::Ptr t(new Track<float>());
  t->set_extends(0, 100);
  t->set_resolution(1);
  t->set_trackname("test");
  t->set_subtrackname("test1");
  for (int i = 0; i < 100; ++i) {
    t->set(i, (float)rand());
  }
  tio_.WriteSubTrack<float>(*t);

  Track<float>::Ptr t2(new Track<float>());
  t->set_extends(0, 100000000);
  t->set_resolution(10);
  t->set_abs_extends(0, 1000000000);
  t->set_trackname("test");
  t->set_subtrackname("test2");
  for (int i = 0; i < 100000000; ++i) {
    t->set(i, (float)rand());
  }
  tio_.WriteSubTrack<float>(*t);
}

TEST_F(HDFTrackTest, ReadTrackTest) {
  Track<float>::Ptr t(new Track<float>());
  tio_.ReadSubTrack(std::string("test"), std::string("test2"),
                    *t);
}

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
