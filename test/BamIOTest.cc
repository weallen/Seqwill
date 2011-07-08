#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <map>
#include <vector>
#include <string>

#include "gtest/gtest.h"

#include "io/BamIO.h"

using namespace std;
namespace {
class BamIOTest : public ::testing::Test {
  protected:
    BamIOTest() {
      fname_ = "/media/Storage/user/data/bam/1_1_sort.bam";
    }

    virtual ~BamIOTest() {
    }

    virtual void SetUp(){
      b_.Init(fname_);
    }
    virtual void TearDown() {
    }

  std::string fname_;
  BamIO b_;
};



TEST_F(BamIOTest, OpenBamTest) {
  ASSERT_TRUE(b_.is_open());
}

TEST_F(BamIOTest, RefSeqInfo) {
  b_.LoadRefSeqInfo();
  ASSERT_EQ(b_.chrlens().size(),22);
  ASSERT_EQ(b_.chrlens()["chr1"], 197195432);
}

TEST_F(BamIOTest, SingleReadFactoryTest) {
  SingleReadFactory f;
  f.Init(b_.ref_seqs());
  BamTools::BamReader r;
  r.Open(fname_);
  BamTools::BamAlignment a;
  r.GetNextAlignment(a);
  f.AddRead(a);
  ASSERT_EQ(f.num_reads(), 1);
  SingleReadFactory::iterator it = f.begin();
  ASSERT_EQ(it->second.pos, a.Position);
  int count = 0;
  while (r.GetNextAlignment(a)) {
    if ((count++ % 1000000) == 0) {
      std::cout << count << " reads..." << std::endl;
    }
    f.AddRead(a);
  }
  std::cout << f.num_reads() << std::endl;
  std::cout << f.MemUsage() << std::endl;
  std::vector<SingleRead> chr1reads = f.GetReadsAtRefID(HashString("chr1"));
  std::cout << chr1reads.size() << std::endl;
}

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
