#ifndef BAM_H_
#define BAM_H_

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <algorithm>

#include <boost/shared_ptr.hpp>

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/SamSequenceDictionary.h>
#include <api/SamSequence.h>
#include "base/Common.h"
#include "base/StringUtil.h"
#include "base/Log.h"
#include "base/FileUtil.h"
#include "base/StringUtil.h"

#include "common/Track.h"


class ReadTable;
class SingleReadFactory;
struct SingleRead;
class BamIO;
typedef uint64_t RefID;
typedef uint64_t ReadID;

typedef enum Stand {
  kFwd,
  kRev,
  kUnknown
} Strand;


struct SingleRead
{
  ReadID insert_id;
  RefID ref_id;
  int32_t pos;

  RefID partner_ref;
  int32_t partner_pos;

  Strand strand;
  // from cigar string
  int32_t len;
  bool is_first;
};

// -------------------------------------------

class ReadHit
{
 public:

  ReadHit() 
    : has_mate_(false) 
    {}

  ReadHit(SingleRead* left_align)
    : has_mate_(false)
    , first_read_(left_align)
    , second_read_(NULL)

  {}

  ReadHit(SingleRead* left_align,
          SingleRead* right_align)
    : has_mate_(true)
    , first_read_(left_align)
    , second_read_(right_align)
  {}

  virtual ~ReadHit()
  {
    delete first_read_;
    delete second_read_;
  }

  SingleRead* left_align() const { return first_read_; }
  SingleRead* right_align() const { return second_read_; }
  const bool has_mate() const { return has_mate_; }

 private:
  DISALLOW_COPY_AND_ASSIGN(ReadHit)

  bool has_mate_;
  SingleRead* first_read_;
  SingleRead* second_read_;
};

// -------------------------------------------

// Holds final representation of
// reads with their mates
// Want to be able to iterate over it, and get reads
// By position and chromosome
class ReadTable
{
public:
  typedef boost::shared_ptr<ReadTable> Ptr;
  typedef std::map<uint64_t, ReadHit*> ReadMap;
  typedef std::map<std::string, uint64_t> InvChrMap;
  typedef std::map<uint64_t, std::string> ChrMap;

  ReadTable()
    : num_reads_(0)
  {}
  virtual ~ReadTable() {}

  void Init(const std::vector<std::string>& chrs);
  void AddReads(SingleReadFactory* sreads);
  
private:

  void Clear();

  DISALLOW_COPY_AND_ASSIGN(ReadTable)

  int num_reads_;

  // member vars
  ReadMap reads_;
  ChrMap chrs_;
  InvChrMap inv_chrs_;
};

// ------------------------------------------

// Convert raw BamAlignments into single reads
// Store the single reads by
class SingleReadFactory
{
public:
  typedef boost::shared_ptr<SingleReadFactory> Ptr;
  typedef std::map<ReadID, SingleRead> ReadMap;

  // sorted by of reads organized by chromosome

  typedef ReadMap::iterator iterator;

  SingleReadFactory()
    : num_reads_(0)
  {}

  virtual ~SingleReadFactory() {}

  void Init(const std::vector<std::string>& refseqnames);

  void AddRead(const BamTools::BamAlignment& read);

  std::vector<SingleRead> GetReadsAtPos(RefID refid, int32_t pos) ;
  std::vector<ReadID> GetReadIDsAtPos(RefID refid, int32_t pos) ;
  std::vector<SingleRead> GetReadsAtRefID(RefID refid);
  SingleRead GetReadByID(ReadID id);

  size_t MemUsage() const;

  // Iterator shit
  iterator begin() { return reads_.begin(); }
  iterator end() { return reads_.end(); }

  int num_reads() const { return num_reads_; }

private:
  DISALLOW_COPY_AND_ASSIGN(SingleReadFactory)

  int num_reads_;
  std::map<int, RefID> refseq_ids_;
  ReadMap reads_;
};

// -------------------------------------------

class BamIO
{
 public:
  typedef boost::shared_ptr<BamIO> Ptr;

  BamIO()  {}
  BamIO(const std::string& fname) 
    { Init(fname); }

  virtual ~BamIO() {
    Close();
  }

  void Init(const std::string& fname);

  SingleReadFactory* LoadAllSingleReads();
  SingleReadFactory* LoadChrSingleReads(const std::string& chr);

  bool is_open() const { return isopen_; }

  std::map<std::string, int> chrlens() const { return chrlens_; }
  std::vector<std::string> ref_seqs() const { return refseqs_; }

 private:
  DISALLOW_COPY_AND_ASSIGN(BamIO)
  void LoadRefSeqInfo();

  bool Open(const std::string& fname);
  void Close();

  bool isopen_;
  BamTools::BamReader reader_;  
  std::map<std::string, int> chrlens_;
  std::vector<std::string> refseqs_;
};

#endif
