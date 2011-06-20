
#include "base/CommandLineParser.h"

#include "base/FileParser.h"
#include "base/StringUtil.h"
#include "analysis/DNAVector.h"


int main(int argc, char** argv)
{

  commandArg<string> aStringCmmd("-fasta_in","input fasta sequence");
  commandArg<string> bStringCmmd("-gff","gff file");
  commandArg<string> cStringCmmd("-fasta_out","output fasta sequence");
  commandLineParser P(argc,argv);
  P.SetDescription("Grabs sequence from a gff file.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(cStringCmmd);
  P.parse();

  string fastbIn = P.GetStringValueFor(aStringCmmd);
  string gffFile = P.GetStringValueFor(bStringCmmd);
  string fastbOut = P.GetStringValueFor(cStringCmmd);

  FlatFileParser fParser(gffFile);


  // fastb file
  vecDNAVector seqs;
  seqs.Read(fastbIn);

  vecDNAVector seqsOut;

  fParser.ParseLine();
  while ( ! fParser.IsEndOfFile() )
  {
//     if ( ! fParser.ParseLine() )
//       break;

    string chr1(""), name1("");
    int start1(0), stop1(0), len1(0), orient1(0);

    string chr2(""), name2("");
    int start2(0), stop2(0), len2(0), orient2(0);

    chr1 = fParser.AsString(0);
    start1 = fParser.AsInt(3);
    stop1 = fParser.AsInt(4);
    len1 = fParser.AsInt(5);
    orient1 = fParser.AsInt(6);
    name1 = fParser.AsString(8);
    
    while ( !fParser.IsEndOfFile() && fParser.AsString(8) == name1 )
    {

      chr2 = fParser.AsString(0);
      start2 = fParser.AsInt(3);
      stop2 = fParser.AsInt(4);
      len2 = fParser.AsInt(5);
      orient2 = fParser.AsInt(6);
      name2 = fParser.AsString(8);
      
      fParser.ParseLine();
    }


    int fastbIndex(-1);
    string chrSub = After(chr1, "chr");
    if ( chrSub == "X" )
      fastbIndex = 23;
    else if ( chrSub == "Y" )
      fastbIndex = 24;
    else
      fastbIndex = atoi(chrSub.c_str());

    cout << chr1 <<" "<< start1 <<" "<< stop1 << " " << len1 <<" "<< fastbIndex << endl;
    cout << chr2 <<" "<< start2 <<" "<< stop2 << " " << len2 <<" "<< fastbIndex << endl;

    DNAVector sub1,sub2;
    sub1.SetToSubOf(seqs[fastbIndex],start1,len1);
    sub2.SetToSubOf(seqs[fastbIndex],start2,len2);
    
    if ( orient1 < 0 )
      sub1.ReverseComplement();
    
    seqsOut.push_back(sub1);
    
    if ( orient2 < 0 )
      sub2.ReverseComplement();
    
    seqsOut.push_back(sub2);
    

  }
  
  seqsOut.Write(fastbOut);
  
}
