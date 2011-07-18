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
class BamIOTest : public ::testing::Test {
  protected:
    BamIOTest() {
    }

    virtual ~BamIOTest() {
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
  

};
TEST_F(BamIOTest, ChromosomeSaveTest) {
  ASSERT_TRUE(SaveChrFromFASTA("/tmp/out.h5", "/media/Storage/user/data/genomedata/fasta/chr1.fa", "mm9"));
}

TEST_F(BamIOTest, ChromosomeLoadTest) {
  Chromosome::Ptr chr(new Chromosome);
  ASSERT_TRUE(LoadChr("/tmp/out.h5", "mm9", "chr1", chr));
  ASSERT_FALSE(LoadChr("/tmp/out.h5", "mm9", "chr2", chr));
  DNASequence::Ptr data = chr->data();
  ASSERT_EQ(data->GetNuc(1), 'N');

}


}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
