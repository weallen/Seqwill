#include <string>

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "analysis/HMMDistance.h"
#include "analysis/DNAVector.h"


////////// Input file format //////////
// Header: name1  name2 ...
// Data:   chr  pos  genotype1  genotype2 ...
// Genotypes: single or double nucleotide
// Missing genotype:  'N'
///////////////////////////////////////



double SimpleDist(char v1, char v2) 
{
  if(v1=='N' || v2=='N' || v1=='X' || v2=='X') return -1.;
  else return DNA_Diff(v1,v2);
}


class SeqPart
{
public:
  SeqPart(const string & name) {
    m_name = name;
  }

  void Set(const string & seq) {
    m_seq = seq;
  } 

  char Get(int i) const {
    if (m_seq == "")
      return '?';
    const char * p = m_seq.c_str();
    return p[i];
  }

  void Reset() {
    m_seq = "";
  }

  const string & Name() const {return m_name;}
private:
  string m_name;
  string m_seq;
};







int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-i","input");
  commandArg<string> bStringCmmd("-o","output");
  commandLineParser P(argc,argv);
  P.SetDescription("Converts data into a HMM-digestable file of features");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
 
  
  P.parse();

  string in = P.GetStringValueFor(aStringCmmd);
  string out = P.GetStringValueFor(bStringCmmd);
  

  HMMFeatureVector f;

  FlatFileParser parser;
  parser.Open(in);


  // get names from first line
  parser.ParseLine();
  svec<SeqPart> seqs;
  int n = parser.GetItemCount();
  if(n==0) {
    cerr << "Wrong format: list names in first line" << endl;
    exit(-1);
  }
  for(int c=0; c<n; c++) {
    f.SetName(c, parser.AsString(c));
    seqs.push_back(SeqPart(parser.AsString(c)));
  }

  int k = 0;
  
  HMMFeature feature;
  
  int i, j, l;
    
  feature.resize(n * n);

  char * line = new char[n];

  string chr;
  int pos = -1;
  string genotype;
  char code;

  while (parser.ParseLine()) {
    
    if (parser.GetItemCount() == 0)
      continue;

    chr = parser.AsString(0);
    pos = parser.AsInt(1);

    for (i=0; i<n; i++) {
      genotype = parser.AsString(i+2);
      
      if(genotype == "-" || genotype == "--"  || genotype == "het" || genotype == "?" || genotype == "N" || genotype == "X")
	code = 'N';
      else  
        code = GetAmbiguous(genotype);
      
      line[i] = code;
    }
          
    if (k >= f.isize())
      f.resize(k + 5000000);
      
    f.SetCompressed(k, line, n, chr, pos);
    k++;
  
      
    for (i=0; i<n; i++)
      seqs[i].Reset();
      
    continue;
   
  }


  f.resize(k);


  f.Write(out);

  delete [] line;

  return 0;

}
  
