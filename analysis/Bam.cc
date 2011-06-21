#include "analysis/Bam.h"

void BamLoader::Close()
{
  reader_.Close();
}

void BamLoader::Init(Genome* g, const std::string& fname)
{
  genome_ = g;
  extend_ = 0;
  if (reader_.IsOpen()) 
    reader_.Close();
  reader_.Open(fname);
  if (!reader_.LocateIndex()) {
    reader_.CreateIndex();
    std::cout << "Creating index for " << fname << "..." << std::endl;
  }
  InitChrs();
}

// Assumes for now that paired reads are on same chromosome
void BamLoader::ReadBamFile()
{
  if (!reader_.IsOpen())
    std::cerr << "Must init BamReader before read file." << std::endl;
  BamAlignment align;
  while (reader_.GetNextAlignment(align)) {
    if (align.IsPaired()) {
      PairedRead(align);
    } else {
      UnpairedRead(align);
    }    
  }
}

void BamLoader::InitChrs()
{
  std::string name;
  int len;
  std::vector<Chromosome*> chrs = genome_->GetChromosomes();
  for (std::iterator<Chromosome*> i = chrs.begin();
       i != chrs.end(); ++i) {
    name = i->GetName();
    len = i->GetLength();
    ublas::zero_vector<float> v(len);
    chrs_[name] = v;
  }
}

void BamLoader::PairedRead(const BamAlignment& read)
{

}

void BamLoader::UnpairedRead(const BamAlignment& read)
{
  int pos;
  
}
