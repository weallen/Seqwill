#ifndef GFF
#define GFF


#include "base/StringUtil.h"
#include "base/FileParser.h"

// note that this class differs from the standard definition of a gff file in that the 
// end base is one past the last base to be included in the feature
class GFFelement 
{
public:
GFFelement() {}
GFFelement(string &name,
	   string &source,
	   string &feature,
	   int start,
	   int end,
	   double score,
	   string &strand,
	   string &frame,
	   string &group)
  : m_name(name), 
    m_source(source),
    m_feature(feature),
    m_start(start),
    m_end(end),
    m_score(score),
    m_strand(strand),
    m_frame(frame),
    m_group(group) {}


int Size() { return m_end-m_start; }

string GetName() const { return m_name; }
int GetStart() const { return m_start; }
int GetEnd() const  { return m_end; }
string GetStrand() const  { return m_strand; }
double GetScore() const { return m_score; }
string GetFeature() const  { return m_feature; }
string GetGroup() const { return m_group; }
string GetSource() const {return m_source; }
string GetFrame() const { return m_frame; }

void SetName(string name) { m_name =name;}
void SetStart(int start) {m_start=start; }
void SetEnd(int end) {m_end=end; }
void SetScore(double s) {m_score=s;}
void SetFeature(string s) {m_feature=s;}
void SetOrientation(string orient) {m_strand=orient;}
void SetGroup(string g) { m_group=g;}
void SetSource(string g) {m_source=g;}
void SetFrame(string g) { m_frame=g; }

bool Empty() { return m_name.empty(); }

bool friend operator < (const GFFelement &l, const GFFelement &r)
{
  if (l.m_name < r.m_name)
    return true;
  if (l.m_name > r.m_name)
    return false;
  if (l.m_start<r.m_start)
    return true;
  if (l.m_start>r.m_start)
    return false;
  if (l.m_end<r.m_end)
    return true;
  if (l.m_end>r.m_end)
    return false;
  if (l.m_score<r.m_score)
    return true;

  return false;
}

bool friend operator == (const GFFelement &l, const GFFelement &r)
{
  return ( l.m_name == r.m_name &&
	   l.m_start==r.m_start &&
	   l.m_end==r.m_end &&
	   l.m_score==r.m_score );

}


void Print(ostream &ostrm)
{
  ostrm << m_name <<"\t"
	<< m_source <<"\t"
	<< m_feature <<"\t"
	<< m_start <<"\t"
	<< m_end <<"\t"
	<< m_score <<"\t"
	<< m_strand <<"\t"
	<< m_frame <<"\t"
	<< m_group << endl;
}

friend ostream &operator << (ostream &ostrm, const GFFelement &e)
{
  ostrm << e.m_name <<"\t"
	<< e.m_source <<"\t"
	<< e.m_feature <<"\t"
	<< e.m_start <<"\t"
	<< e.m_end <<"\t"
	<< e.m_score <<"\t"
	<< e.m_strand <<"\t"
	<< e.m_frame <<"\t"
	<< e.m_group;

  return ostrm;
}

bool Overlaps(GFFelement &e)
{
  if (m_name != e.GetName())
    return false;

  int beg = std::max(m_start,e.GetStart());
  int end = std::min(m_end,e.GetEnd());
  
  if (end-beg<0)
    return false;

  return true;
}


private:
string m_name, m_source, m_feature;
int m_start,m_end;
double m_score;
string m_strand, m_frame, m_group;

};


struct orderByGroup : public binary_function<const GFFelement, const GFFelement, bool>
{
bool operator() (const GFFelement &l, const GFFelement &r) const
{
  return (l.GetGroup()<r.GetGroup());
}

};

struct orderByName : public binary_function<const GFFelement, const GFFelement, bool>
{
bool operator() (const GFFelement &l, const GFFelement &r) const
{
  return (l.GetName()<r.GetName());
}

};


inline void GFFWindow(vector<GFFelement> &in, int w, vector<int> &counts, vector<int> &starts)
{
  counts.clear();
  starts.clear();
  if (in.empty())
    return;

  int num_windows = (in[in.size()-1].GetEnd()-in[0].GetStart())/w + 2;

  counts.resize(num_windows), starts.resize(num_windows);

  int this_start = std::max(in[0].GetStart()-w,0);
  int this_end = this_start+w;
  int this_window(0);


    int indx=0;

  while (this_window<num_windows)
  {
    starts[this_window]=this_start;

    while (indx<(int) in.size() && in[indx].GetEnd()<this_end)
    {
      ++counts[this_window];
      ++indx;
    }
    
    ++this_window;
    this_start += w;
    this_end += w;
    
  }
}



void LoadGFFfile(string &filename,vector<GFFelement> &elements)
{
  elements.clear();
  
  FlatFileParser fp;
  
  fp.Open(filename);
  while (fp.ParseLine())
  {
    string name=fp.AsString(0);
    string sour=fp.AsString(1);
    string fea=fp.AsString(2);
    int start = fp.AsInt(3);
    int end = fp.AsInt(4);
    double sc= fp.AsFloat(5);
    string strd = fp.AsString(6);
    string frm=fp.AsString(7);
    string gp = fp.AsString(8);

    elements.push_back(GFFelement(name,sour,fea,start,end,sc,strd,frm,gp));

  }
}



void LoadGencode(string &filename,vector<GFFelement> &elements)
{

  elements.clear();

  FlatFileParser fp;
  fp.Open(filename);
  string typetx("transcript"), typeex("coding"), type5putr("utr5"), type3putr("utr3");
  string typeup150("up200"), typeup2k("up2000"), typefi("first_intron"), typein("intron"), typenc("non-coding"), typeing("intergenic");

  string delim(","), bit("#");
  while (fp.ParseLine())
  {
    string b(fp.AsString(0));
    if (Contains(b,bit))
      continue;

    int bin = fp.AsInt(0);
   
    string orient("1");
    string gene_name = fp.AsString(1);
    string chr = fp.AsString(2);
    string strand = fp.AsString(3); // + or -
    if (strand=="-")
      orient = "-1";

    double score = fp.AsFloat(11);
    
    int txStart = fp.AsInt(4);
    int txEnd = fp.AsInt(5);
    int cdsStart=fp.AsInt(6);
    int cdsEnd = fp.AsInt(7);
    int exonCount = fp.AsInt(8);


    string exonStarts = fp.AsString(9);
    string exonEnds = fp.AsString(10);
    string gene_name2 = fp.AsString(12);

    StringParser sp, ep;
    sp.SetLine(exonStarts,delim);
    int exon_count=sp.GetItemCount();
    if (exon_count != exonCount) {
      cout << "Yikeie "<< exon_count <<" "<< exonCount << " "<< gene_name << endl;
      exit(-1);
    }
    ep.SetLine(exonEnds,delim);
    int exon_counte=ep.GetItemCount();
    if (exon_counte != exonCount) {
      cout << "Yikeieeee"<<endl;
      exit(-1);
    }

    bool isrc(false);
    if (orient == "-1")
      isrc=true;

    bool isnc(false);
    
    string tmptype = typetx;

    if (cdsStart==txEnd) {
      isnc=true;
      tmptype=typenc;
    }

    // transcript
    elements.push_back(GFFelement(chr,
				  tmptype,
				  tmptype,
				  txStart,
				  txEnd,
				  score,
				  orient,
				  gene_name2,
				  gene_name));

    if (isnc)
      continue;

    
    // 5' utr
    if (!isrc && txStart != cdsStart)
      elements.push_back(GFFelement(chr,
				    type5putr,
				    type5putr,
				    txStart,
				    cdsStart,
				    score,
				    orient,
				    gene_name2,
				    gene_name));
    else if (isrc && cdsEnd != txEnd)
      elements.push_back(GFFelement(chr,
				    type5putr,
				    type5putr,
				    cdsEnd,
				    txEnd,
				    score,
				    orient,
				    gene_name2,
				    gene_name));
    
    // 3' utr
    if (!isrc && cdsEnd != txEnd)
      elements.push_back(GFFelement(chr,
				    type3putr,
				    type3putr,
				    cdsEnd,
				    txEnd,
				    score,
				    orient,
				    gene_name2,
				    gene_name));
    else if (isrc && txStart != cdsStart)
      elements.push_back(GFFelement(chr,
				    type3putr,
				    type3putr,
				    txStart,
				    cdsStart,
				    score,
				    orient,
				    gene_name2,
				    gene_name));
  
   
    
    
    // up150
    int up150start = std::max(0,txStart-200);
    int up150end = txStart;
    if (isrc) {
      up150start = txEnd;
      up150end = up150start+200;
    }
    elements.push_back(GFFelement(chr,
				  typeup150,
				  typeup150,
				  up150start,
				  up150end,
				  score,
				  orient,
				  gene_name2,
				  gene_name));
    
    // up2000
    int up2kstart = std::max(0,txStart-2000);
    int up2kend = txStart;
    if (isrc) {
      up2kstart = txEnd;
      up2kend = up2kstart+2000;
    }
    elements.push_back(GFFelement(chr,
				  typeup2k,
				  typeup2k,
				  up2kstart,
				  up2kend,
				  score,
				  orient,
				  gene_name2,
				  gene_name));
    
    
    

    // get the exons and introns
    for (int k=0; k<exonCount; ++k)
    {
      int this_start = std::max(cdsStart,sp.AsInt(k));
      int this_end = std::min(cdsEnd,ep.AsInt(k));

      bool need_exon(true);
      if (k==0 && this_end<cdsStart)
	need_exon=false;

      if (k==exonCount-1 && this_start>cdsEnd)
	need_exon=false;

      if (need_exon) 
      {
	elements.push_back(GFFelement(chr,
				      typeex,
				      typeex,
				      this_start,
				      this_end,
				      score,
				      orient,
				      gene_name2,
				      gene_name));
      }

      // introns
      if (k<exonCount-1)  
      {
	this_start = ep.AsInt(k);
	this_end = sp.AsInt(k+1);
	string tmptype = typein;
	if ( (!isrc && k==0) || (isrc && k==exonCount-2))
	  tmptype = typefi;

	elements.push_back(GFFelement(chr,
				      tmptype,
				      tmptype,
				      this_start,
				      this_end,
				      score,
				      orient,
				      gene_name2,
				      gene_name));
      }
    }

  }
}



#endif
