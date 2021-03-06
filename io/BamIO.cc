#include "io/BamIO.h"

void
ReadTable::Init(const std::vector<std::string>& chrnames)
{
  uint32_t currid;

  for (size_t i = 0; i < chrnames.size(); ++i) {
    currid = HashString(chrnames[i].c_str());
    chrs_[currid] = chrnames[i];
    inv_chrs_[chrnames[i]] = currid;
  }
}

// Goes through SingleReads, finds the mate for each,
// And combines them into a single ReadHit
void
ReadTable::AddReads(SingleReadFactory* sreads)
{
  std::map<uint32_t, SingleRead> first_reads;
  std::map<uint32_t, SingleRead> second_reads;
  SingleRead currs;

  for (SingleReadFactory::iterator i = sreads->begin();
       i != sreads->end(); ++i) {
    currs = static_cast<SingleRead>(i->second);
  }
}

void
ReadTable::Clear()
{
  for (ReadMap::iterator i = reads_.begin();
       i != reads_.end(); ++i) {
    delete i->second;
  }
}

// -------------------------------------------

void
SingleReadFactory::Init(const std::vector<std::string>& refseqnames)
{
  for (size_t i = 0; i < refseqnames.size(); ++i) {
    refseq_ids_[i] = HashString(refseqnames[i].c_str());
  }
}

std::vector<SingleRead>
SingleReadFactory::GetReadsAtPos(RefID refid, int32_t pos)
{
  std::vector<SingleRead> reads;
  for (ReadMap::iterator it = reads_.begin();
       it != reads_.end(); ++it) {
    if (it->second.pos == pos && it->second.ref_id == refid)
      reads.push_back(it->second);
  }
  return reads;
}

std::vector<ReadID>
SingleReadFactory::GetReadIDsAtPos(RefID refid, int32_t pos)
{
  std::vector<ReadID> readids;
  for (ReadMap::iterator it = reads_.begin();
       it != reads_.end(); ++it) {
    if (it->second.pos == pos && it->second.ref_id == refid)
      readids.push_back(it->first);
  }
  return readids;
}


SingleRead
SingleReadFactory::GetReadByID(ReadID id)
{
  return reads_[id];
}

size_t
SingleReadFactory::MemUsage() const
{
  return sizeof(*this);
}


std::vector<SingleRead> 
SingleReadFactory::GetReadsAtRefID(RefID id)
{
  std::vector<SingleRead> reads;
  for (ReadMap::iterator it = reads_.begin();
       it != reads_.end(); ++it) {
    if (it->second.ref_id == id)
      reads.push_back(it->second);
  }
  return reads;
}


void
SingleReadFactory::AddRead(const BamTools::BamAlignment& read)
{
  if (read.IsMapped()) {
    SingleRead r;
    r.insert_id = HashString(read.Name.c_str());
    r.ref_id = refseq_ids_[read.RefID];
    r.pos = read.Position;
    r.len = read.Length;
    r.is_first = read.IsFirstMate();
    if (read.IsReverseStrand()) {
      r.strand = kRev;
    } else {
      r.strand = kFwd;
    }
    if (read.IsPaired() && read.IsMateMapped()) {
      r.partner_ref = refseq_ids_[read.MateRefID];
      r.partner_pos = read.MatePosition;
    } else {
      r.partner_pos = -1;
      r.partner_ref = -1;
    }
    reads_[r.insert_id] = r;
    num_reads_++;
    //AddReadToRefMap(r);
  }
}

// -------------------------------------------

bool BamIO::Open(const std::string& fname)
{
  if (!FileExists(fname)) {
    ERRORLOG("File doesn't exist " + fname);
    return false;
  }
  reader_->Open(fname);
  isopen_ = true;
  return true;
}

void BamIO::Close()
{
  if (isopen_) {
    isopen_ = false;
    reader_->Close();
  }
}

void BamIO::Init(const std::string& fname)
{
  if (!reader_->IsOpen())
    Open(fname);
  if (!reader_->LocateIndex()) {
    reader_->CreateIndex();
    DEBUGLOG("Creating index for " + fname + "...");
  }
  LoadRefSeqInfo();
}

// Assumes for now that paired reads are on same chromosome
SingleReadFactory*
BamIO::LoadAllSingleReads()
{
  SingleReadFactory* factory = new SingleReadFactory;
  factory->Init(refseqs_);
  BamTools::BamAlignment align;
  for (std::vector<std::string>::const_iterator it = refseqs_.begin();
       it != refseqs_.end(); ++it) {
    int refid = reader_->GetReferenceID(*it);
    reader_->SetRegion(refid, 0, refid, chrlens_[*it]);
    while (reader_->GetNextAlignment(align)) {
      factory->AddRead(align);
    }
    reader_->Rewind();    
  }
  return factory;
}

SingleReadFactory* 
BamIO::LoadChrSingleReads(const std::string& chr)
{
  SingleReadFactory* factory = new SingleReadFactory;
  BamTools::BamAlignment align;
  factory->Init(refseqs_);
  int refid = reader_->GetReferenceID(chr);
  reader_->SetRegion(refid, 0, refid, chrlens_[chr]);
  while (reader_->GetNextAlignment(align)) {
    factory->AddRead(align);
  }
  reader_->Rewind();
  return factory;
}

void 
BamIO::LoadRefSeqInfo()
{
  chrlens_.clear();
  refseqs_.clear();
  BamTools::SamHeader header = reader_->GetHeader();
  if (header.HasSequences()) {
    BamTools::SamSequenceDictionary seqs = header.Sequences;
    BamTools::SamSequence currseq;
    for (BamTools::SamSequenceIterator i = seqs.Begin();
         i != seqs.End(); ++i) {
      currseq = *i;
      if (currseq.HasName() && currseq.HasLength()) {
          chrlens_[currseq.Name] = StringToInt(currseq.Length);
          refseqs_.push_back(currseq.Name);
      }
    }   
  }
}
