#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <string>

#include "gtest/gtest.h"
#include "analysis/Chromosome.h"
#include "base/SVector.h"
#include "analysis/Genome.h"

using namespace std;
namespace {
class ChromosomeTest : public ::testing::Test {
  protected:
    ChromosomeTest() 
      : m_genome("testdata")
    {
    }

    virtual ~ChromosomeTest() {
      system("rm -rf testdata");
    }

    virtual void SetUp() {
      m_genome.Open();
    }
    virtual void TearDown() {
      m_genome.Close();
    }
  
    string m_hdf5_dir;
    string m_fasta_dir;
    Genome m_genome;
};

TEST_F(ChromosomeTest, OpenCloseTest) {
  Chromosome c("testchr", "testdata");
  c.Open();
  ASSERT_TRUE(c.IsOpen());
  ASSERT_EQ(c.GetName(), "testchr");
  ASSERT_EQ(c.Close(), 0);
  ASSERT_FALSE(c.IsOpen());
}

TEST_F(ChromosomeTest, ReadWriteSeqTest) {
  Chromosome c("testchr", "testdata");
  string s("GATTATACCCTTATATATCTATATATACTATA");
  ASSERT_GE(c.WriteSeq(s), 1);
  string s2;
  ASSERT_GE(c.ReadSeq(&s), 1);
  ASSERT_EQ(s2, s);
}

TEST_F(ChromosomeTest, SeqLenTest) {
  Chromosome c("testchr", "testdata");
  string s("GATTATACCCTTATATATCTATATATACTATA");
  ASSERT_EQ(c.GetLength(), 0);
  c.WriteSeq(s);
  ASSERT_EQ(c.GetLength(), 32);
  ASSERT_EQ(c.GetLength(), 32);
}

TEST_F(ChromosomeTest, ReadWriteTrackTest) {
  ublas::vector<float> v(100);
  for (int i=0; i<100; i++) {
    v(i) = 1;
  }
  string chrname("chr10");
  Chromosome* chr;
  m_genome.GetChromosome(chrname, chr);
  chr->WriteTrack("track1", v);
  //ASSERT_EQ(chr->TrackNames()[0], "track1");
}

TEST_F(ChromosomeTest, DeleteTrackTest) {
  
}

} // namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
