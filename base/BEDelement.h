
#ifndef BEDELEMENT_H
#define BEDELEMENT_H

#include <string>
#include <sstream>
#include "base/StringUtil.h"
#include "base/FileParser.h"
#include <iostream>
#include <map>
#include <vector>


class BEDelement
{
public:
BEDelement () {}

BEDelement(string chr, int start, int end, string name)
  : m_chr(chr), m_start(start), m_end(end), m_name(name) {}


string GetChromosome() const {return m_chr; }
string GetName() const {return m_name;}
int GetStart() const {return m_start;}
int GetEnd() const {return m_end;}

void SetEnd(int e) {m_end=e;}
void SetStart(int s) {m_start=s;}

int Size() {return m_end-m_start;}

bool Subsumes(BEDelement &e)
{
  if (m_chr != e.GetChromosome())
    return false;

  return (m_start<=e.GetStart() && m_end>=e.GetEnd());
}

bool IsSubsumedBy(BEDelement &e)
{
  if (m_chr != e.GetChromosome())
    return false;

  return (e.GetStart()<=m_start && e.GetEnd()>=m_end);
}

bool Overlaps(BEDelement &e)
{
  if (m_chr != e.GetChromosome())
    return false;

  int beg = std::max(m_start,e.GetStart());
  int end = std::min(m_end,e.GetEnd());
  
  if (end-beg<0)
    return false;

  return true;
}

friend ostream &operator << (ostream &ostrm, const BEDelement &e)
{
  ostrm << e.m_chr <<"\t"
	<< e.m_start<<"\t"
	<< e.m_end<<"\t"
	<< e.m_name;

  return ostrm;
} 

bool friend operator < (const BEDelement &l, const BEDelement &r)
{
  if (l.m_chr < r.m_chr)
    return true;
  if (l.m_chr > r.m_chr)
    return false;
  if (l.m_start<r.m_start)
    return true;
  if (l.m_start>r.m_start)
    return false;
  if (l.m_end<r.m_end)
    return true;
  if (l.m_end>r.m_end)
    return false;
  if (l.m_name<r.m_name)
    return true;

  return false;
}

bool friend operator == (const BEDelement &l, const BEDelement &r)
{
  if (l.m_chr != r.m_chr)
    return false;
  if (l.m_start!=r.m_start)
    return false;
  if (l.m_end!=r.m_end)
    return false;
  if (l.m_name != r.m_name)
    return false;

  return true;
}



private:
string m_chr;
int m_start,m_end;
string m_name;

};


void LoadBEDfile(string &filename, vector<BEDelement> &elements);





// int Size(vector<BEDelement> &b,string &chr)
// {
//   int indx(0);
//   int num(b.size());
//   while (indx<num && b[indx].GetChromosome() != chr)
//     ++indx;

//   int tot(0);
//   while (indx<num && b[indx].GetChromosome() == chr)
//   {
//     tot += b[indx].Size();

//     ++indx;
//   }
  
//   return tot;

// }

// int UniqueSize(vector<BEDelement> &b, string &chr)
// {
//   int indx(0);
//   int num(b.size());
//   while (indx<num && indx<num && b[indx].GetChromosome() != chr)
//     ++indx;

//   int min = b[indx].GetStart();
//   int max = b[indx].GetStart();

//   while (indx<num && b[indx].GetChromosome() == chr)
//   {
//     if (b[indx].GetStart()<min)
//       min=b[indx].GetStart();
//     if (b[indx].GetEnd()> max)
//       max = b[indx].GetEnd();

//     ++indx;
//   }

//   vector<int> range(max-min,0);
//   indx=0;
//   while (indx<num && b[indx].GetChromosome() != chr)
//     ++indx;

//   while (b[indx].GetChromosome() == chr)
//   {
//     for (int i=b[indx].GetStart(); 
//          i<b[indx].GetEnd(); ++i)
//       range[i-min]=1;

//     ++indx;
//   }

//   int tot(0);
//   for (int i=0; i<(int) range.size(); ++i)
//     tot += range[i];

//   return tot;
// }


#endif
