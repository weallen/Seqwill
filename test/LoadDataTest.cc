#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <vector>
#include <string>

#include "gtest/gtest.h"

#include "base/Types.h"
#include "analysis/Genome.h"

using namespace std;
namespace {
class LoadDataTest : public ::testing::Test {
  protected:
    LoadDataTest() {
      filename_ = "/Users/admin/src/Seqwill/testdata/moe_d3a_hmc_50_raw.wig";
      outname_ = "/Users/admin/Documents/test.trk";
      genome_info_name_ = "/Users/admin/Documents/chromosomes.trk";
    }

    virtual ~LoadDataTest() {
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
  std::string genome_info_name_;
  std::string filename_;
  std::string outname_;
};
      
TEST_F(LoadDataTest, LoadGenomeInfoTest) {
  GenomeInfo g;
  LoadGenomeInfoFromChr(genome_info_name_, std::string("mm9"), &g);
  std::vector<std::string> c = g.chr_names();
  ASSERT_EQ(g.chr_size(std::string("chr1")), 197195432);
}
 

TEST_F(LoadDataTest, WIGToTrackTest) {
   Track<float> t;
   GenomeData d;
   GenomeInfo g;
   LoadGenomeInfoFromChr(genome_info_name_, std::string("mm9"), &g);
   d.Init(outname_, g);
   d.SaveTrackFromWIG(filename_, std::string("moe_d3a_hmc_raw"), 50);
}

TEST_F(LoadDataTest, LoadFromWIGTest) {

  WIGParser p;
  p.Open(filename_);
  std::vector<WIGLine> lines;
  for (int i = 0; i < 1000000; ++i) {
    lines.push_back(p.NextLine());
  }
  ASSERT_EQ(p.num_lines(), 53435119);
  ASSERT_EQ(lines[0].val, 0.0);
  ASSERT_EQ(p.curr_state().chr, kChr1);
  
  for (int i = 0; i < 5000000; ++i) {
    p.NextLine();
  }
  ASSERT_EQ(p.curr_state().chr, kChr2);
  
  //ASSERT_EQ(0.778794, lines[100000].val);
}

}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
