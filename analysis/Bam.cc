#include "analysis/Bam.h"

void BamLoader::Close()
{
  reader_.Close();
}

void BamLoader::Init(TrackReader* g, const std::string& fname)
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
  BamTools::BamAlignment align;
  while (reader_.GetNextAlignment(align)) {
    if (align.IsPaired()) {
      PairedRead(align);
    } else {
      UnpairedRead(align);
    }    
  }
}

std::map<std::string, int> BamLoader::GetRefSeqInfo() 
{
  std::vector<std::string> refseqs;
  BamTools::SamHeader header = reader_.GetHeader();
  if (header.HasSequences()) {
    BamTools::SamSequenceDictionary seqs = header.Sequences();
    BamTools::SamSequence currseq;
    for (BamTools::SamSequenceIterator i = seqs.Begin();
	 i != seqs.End(); ++i) {
      currseq = *i;
      if (currseq.HasName() && currseq.HasLength())
	refseqs.insert(std::make_pair(currseq.Name(), currseq.Length()));
    }   
  }
  return refseqs;
}

void BamLoader::InitChrs()
{
  std::string name;
  int len;
  std::vector<Chromosome*> chrs = genome_->GetChromosomes();
  for (std::vector<Chromosome*>::iterator i = chrs.begin();
       i != chrs.end(); ++i) {
    name = (*i)->GetName();
    len = (*i)->GetLength();    
    chrs_[name] = v;
  }
}

void BamLoader::PairedRead(const BamTools::BamAlignment& read)
{

}

void BamLoader::UnpairedRead(const BamTools::BamAlignment& read)
{
  int pos;
  
}
