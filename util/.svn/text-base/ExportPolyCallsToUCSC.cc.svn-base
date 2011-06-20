
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <string>
#include <sstream>

struct cv_track
{
cv_track() {}

cv_track(int s,int e, double sc, string c)
  : startpt(s), endpt(e), score(sc), chr(c) {}


bool operator < (const cv_track &c) const
{
  if ( chr < c.chr )
    return true;
  if ( chr > c.chr )
    return false;

  if ( startpt < c.startpt )
    return true;

  return false;

}

int startpt,endpt;
double score;
string chr;
};




struct poly_call
{

poly_call() {}

poly_call(int rid, int rpos, int dcov, int acov, int tcov, int gcov, int ccov)
  : refId(rid), refPos(rpos), readCov(dcov), aCov(acov), tCov(tcov), gCov(gcov), cCov(ccov) {}


bool operator < (const poly_call &c) const
{
  if ( refId < c.refId )
    return true;
  if ( refId > c.refId )
    return false;

  if ( refPos < c.refPos )
    return true;
  
  return false;
}

bool operator == (const poly_call &c) const
{
  if ( refId == c.refId && refPos == c.refPos )
    return true;

  return false;
}

int MajorityCount()
{
  return (std::max(aCov,std::max(tCov,std::max(gCov,cCov))));
}

string MajorityBase()
{
  int maxc = std::max(aCov,std::max(tCov,std::max(gCov,cCov)));
  string mb("");

  if ( maxc == aCov )
    mb+="A";
  if ( maxc == tCov )
    mb+="T";
  if ( maxc == gCov )
    mb+="G";
  if ( maxc == cCov )
    mb+="C";

  return mb;
}


int refId, refPos, readCov, aCov, tCov, gCov, cCov, refMatch;
string refBase;
bool isHet;

};

void FillCalls(string &file, vector<poly_call> &calls)
{
  calls.clear();

  FlatFileParser fp(file);
  string line("");

  int count(0);
  while ( fp.GetLine(line) )
  {
    StringParser sp;
    sp.SetLine(line);

    string::size_type t = line.find("contig");
    if ( t != string::npos )
      continue;

    t=line.find("...");
    if ( t != string::npos )
      continue;

    t=line.find("Version");
    if ( t != string::npos )
      continue;
      
    poly_call p(sp.AsInt(0),
		sp.AsInt(1),
		sp.AsInt(3),
		(sp.AsInt(4)+sp.AsInt(5)),
		(sp.AsInt(6)+sp.AsInt(7)),
		(sp.AsInt(8)+sp.AsInt(9)),
		(sp.AsInt(10)+sp.AsInt(11)));

    p.refMatch=sp.AsInt(12);
    p.isHet=true;
    if ( sp.AsString(13) == "homozygous")
      p.isHet=false;
    p.refBase = sp.AsString(2);


    calls.push_back(p);
  }
}

int main(int argc, char** argv)
{

  commandArg<string> covFile("-f","file of polyCalls output");
  commandArg<string> chr("-c","target chromosome");
  commandArg<int> targId("-t","target id");
  commandArg<int> offset("-o","genomic coord for start of target");
  commandArg<string> trackname("-n","name for track in UCSC");

  commandArg<string> conservFile("-consv","use the UCSC conservation track to filter snps","");
  commandArg<double> minConsvScore("-ms","min conservation score",0.0);


  commandLineParser P(argc,argv);

  P.registerArg(covFile);
  P.registerArg(targId);
  P.registerArg(offset);
  P.registerArg(chr);
  P.registerArg(trackname);
  P.registerArg(conservFile);
  P.registerArg(minConsvScore);

  P.parse();

  string some = P.GetStringValueFor(chr);
  string file  = P.GetStringValueFor(covFile);
  int id = P.GetIntValueFor(targId);
  int start = P.GetIntValueFor(offset);
  string name=P.GetStringValueFor(trackname);

  string cFile = P.GetStringValueFor(conservFile);
  double minScore = P.GetDoubleValueFor(minConsvScore);

  vector<cv_track> cvTrack;
  
  if ( !cFile.empty() )
  {
    FlatFileParser fp(cFile);
    string line("");
    while ( fp.GetLine(line) )
    {
      StringParser sp;
      sp.SetLine(line);
      
      
      if ( line.find("#") != string::npos )
	continue;
      
      if ( sp.AsFloat(5) > minScore )
	continue;
      
      cv_track tmp(sp.AsInt(2), 
		   sp.AsInt(3),
		   sp.AsFloat(5),
		   sp.AsString(1));
      
      cvTrack.push_back(tmp);
    
    }
    sort(cvTrack.begin(),cvTrack.end());
  }

  vector<poly_call> calls;
  FillCalls(file,calls);  
  sort(calls.begin(),calls.end());

//   cout << "%chromosome source SNP start end counts(A,T,G,C,matching reference) coverage snp"<<endl;

  cout << "track name=" << flush;
  cout << name << endl;

  string desc("polyCalls\tSNP");


  for (int i=0; i<(int) calls.size(); ++i )
  {
    if ( calls[i].refId != id )
      continue;

    int snpStart = calls[i].refPos+start;
    int snpEnd = calls[i].refPos+start;

    bool keeper(true);
    
    if ( !cvTrack.empty() )
    {
      desc = "polyCalls\tSNP_in_conserved";

      keeper=false;
      for ( int j=0; j<(int) cvTrack.size(); ++j )
      {
	if ( cvTrack[j].chr != some )
	  continue;

	if ( snpEnd >= cvTrack[j].startpt && snpStart< cvTrack[j].endpt )
	{
	  keeper =true;
	  break;
	}
      }
    }
    
    if ( !keeper )
      continue;
    

    ostringstream ost;
    ost << some <<"\t"<<desc<<"\t";
    ost << snpStart <<"\t"<< snpEnd << "\t";
    ost << "A" << calls[i].aCov 
	<< "_T" << calls[i].tCov
	<< "_G" << calls[i].gCov
	<< "_C" << calls[i].cCov
	<< "_R(" << calls[i].refBase <<")"<< calls[i].refMatch << "_";
    
    string zyg("homozygous");
    if ( calls[i].isHet )
      zyg="heterozygous";
    
    ost << zyg;
    
    cout << ost.str() << "\t.\t"
	 << calls[i].readCov << "\t"
	 << calls[i].MajorityBase()
	 << endl;

  }
}







