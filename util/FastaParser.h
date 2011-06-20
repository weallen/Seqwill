#ifndef _FASTAPARSER_H_
#define _FASTAPARSER_H_

#include "base/FileParser.h"
#include "analysis/DNAVector.h"
#include "math/Functions.h"


class FastaSequence
{
 public:
   FastaSequence() {};
   FastaSequence(string & s, DNAVector & d) {
     m_seq=d;
     m_name=s;
   }

   DNAVector GetSeq() {return m_seq;}
   string GetName() {return m_name;}

   void Print() {
     cout << ">" << m_name << endl;
     m_seq.Print();
   }

   void SetSeq(DNAVector & d) {m_seq=d;}
   void SetName(string & s) {m_name=s;}

 private:
   DNAVector m_seq;
   string m_name;
};




class FastaParser // parse fasta file and make available the names and sequences 
{
 public:
   FastaParser() {};
   
   void Open(string & file) {fileparse.Open(file);}
   int CountSeqs() {return m_numseqs;}
   string GetName(int n) {return m_seqs[n].GetName();}
   DNAVector GetDNAVector(int n) {return m_seqs[n].GetSeq();}
 


   void Parse() {
     m_numseqs = 0;
     FastaSequence currentseq;
     string lastname, currentline, tempseq;
     bool first = 1;
     DNAVector d;
     while(fileparse.GetLine(currentline)) {
       if(currentline[0]=='>') {
	 m_numseqs++;
	 if(!first) {
	   currentseq.SetName(lastname);
	   d.SetFromBases(tempseq);
	   currentseq.SetSeq(d);
	   tempseq.clear();
	   m_seqs.push_back(currentseq);
	 }
	 first = 0;
      	 lastname = currentline;
       } else tempseq.append(currentline);
     }
     d.SetFromBases(tempseq);
     currentseq.SetSeq(d);
     currentseq.SetName(lastname);
     m_seqs.push_back(currentseq);

   }



 
 private:
   FlatFileParser fileparse;
   vec<FastaSequence> m_seqs;
   int m_numseqs;
};




class FastqParser // parse fastq file and make available the names and quals as vector or string
{
 public:
   FastqParser() {};
   
   void Open(string & file) {fileparse.Open(file);}
   
   int CountSeqs() {return m_numseqs;}
   string GetName(int n) {return m_names[n];}
   vector<int> GetQualVec(int n) {return m_quals[n];}
   vector<int> GetQualHistogram(int n) {return m_histograms[n];}
   string GetQualString(int n) {return m_qualstrings[n];}
   int GetLength(int n) {return m_lengths[n];}
   
   void PrintQuals(int n) {
     int l=m_lengths[n];
     for(int i=0; i<l; i++) cout << m_quals[n][i] << " ";
     cout << endl;
   }
   void PrintHistogram(int n) {
     int l=m_histograms[n].size();
     for(int i=0; i<l; i++) cout << m_histograms[n][i] << " ";
     cout << endl;
   }

   double GetQualMean(int n) {
     return Mean(m_quals[n]);
   }

   double GetHistogramMean(int n) {
     return DistMean(m_histograms[n]);
   }

   double GetQualMedian(int n) {
     return Median(m_quals[n]);
   }

   double GetHistogramMedian(int n) {
     return DistMedian(m_histograms[n]);
   }

   double GetQualVariance(int n) {
     return Variance(m_quals[n],Mean(m_quals[n]));
   }

   double GetHistogramVariance(int n) {
     return Variance(m_histograms[n],Mean(m_histograms[n]));
   }

   void Parse() {
       m_numseqs=0;
       string currentline;
       string currentqual;
       while(fileparse.GetLine(currentline)) {
	  if(currentline[0]=='>') {
              m_numseqs++;
	      currentline.erase(0,1);
	      m_names.push_back(currentline);
	      if(!currentqual.empty()) {
	          m_qualstrings.push_back(currentqual);
	          currentqual.clear();
		}  
	  } else {currentqual.append(" "); currentqual.append(currentline);}
	 }
       m_qualstrings.push_back(currentqual);
      
       StringParser getmaxqual, qualstringparse;
       int seqlength,maxqual,tempqual;
       for(int i=0; i<m_numseqs; i++) {
	 getmaxqual.SetLine(m_qualstrings[i]);
	 seqlength=getmaxqual.GetItemCount();
	 m_lengths.push_back(seqlength);
	 maxqual=0;
	 for(int j=0; j<seqlength; j++) {
	   tempqual=getmaxqual.AsInt(j);
	   if(tempqual>maxqual) maxqual=tempqual;
	 }
	 vector<int> histogram(maxqual+1);
	 for(int k=0; k<maxqual+1; k++) histogram[k]=0;
	 vector<int> quals(seqlength);
	 qualstringparse.SetLine(m_qualstrings[i]);
	 for(int l=0; l<seqlength; l++) {
	   tempqual=qualstringparse.AsInt(l);
	   histogram[tempqual]++;
	   quals[l]=tempqual;
	 }
       	 m_quals.push_back(quals);
	 m_histograms.push_back(histogram);
       }
	 
     }
 
 private:
   FlatFileParser fileparse;
   vec<string> m_names;
   vec<string> m_qualstrings;
   vec<int> m_lengths;
   vector< vector<int> > m_quals;
   vector< vector<int> > m_histograms;
   int m_numseqs;
};




#endif 
