#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <string>

#include "gtest/gtest.h"
#include "analysis/Genome.h"
#include "base/SVector.h"

using namespace std;
namespace {
class GenomeTest : public ::testing::Test {
  protected:
    GenomeTest() {  
      m_hdf5_dir = "data"; 
      m_fasta_dir = "/media/Storage/user/data/genomedata/fasta/";
      m_genome = Genome("data");
    }

    virtual ~GenomeTest() {
      system("rm -rf data");
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
  
    string m_hdf5_dir;
    string m_fasta_dir;
    Genome m_genome;
};
TEST_F(GenomeTest, OpenCloseTest) {
  ASSERT_FALSE(m_genome.IsOpen());
  m_genome.Open();
  ASSERT_TRUE(m_genome.IsOpen());
  m_genome.Close();
  ASSERT_FALSE(m_genome.IsOpen());
}

TEST_F(GenomeTest, LoadSeqTest) {
  m_genome.Open();
  if (m_genome.GetAllChromosomeNames().size() > 0) {
    return;
  }
  struct stat filestat;
  DIR* fastaDir = opendir(m_fasta_dir.c_str());
  struct dirent* dirp;
  while ((dirp = readdir(fastaDir))) {
    string fname(dirp->d_name);
    string extension(".fa");
    if (Contains(fname, extension)) {
      cerr << fname << endl;
      m_genome.LoadSeq(m_fasta_dir + "/" + fname);
    }
  }
  svec<string> chrs = m_genome.GetAllChromosomeNames();
  ASSERT_EQ(chrs[0], "chr19");
  m_genome.Close();
}

TEST_F(GenomeTest, GetChrNamesTest) {
  m_genome.Open();
  svec<string> chrnames = m_genome.GetAllChromosomeNames();
  ASSERT_NE(chrnames.size(), 0);
  m_genome.Close();
}

TEST_F(GenomeTest, GetChrTest) {
  m_genome.Open();
  string chrname("chr10");
  ASSERT_TRUE(m_genome.HasChromosome(chrname));
  Chromosome* chr;
  m_genome.GetChromosome(chrname, chr);
  ASSERT_EQ(chr->GetName(), "chr10");
  string chrname2("chrblah");
  ASSERT_FALSE(m_genome.HasChromosome(chrname2));
  m_genome.GetChromosome(chrname2, chr);
  ASSERT_FALSE(chr);
}

TEST_F(GenomeTest, GetTrackNamesTest) {
}

} // namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
