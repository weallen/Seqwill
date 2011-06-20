
#ifndef INTERVAL_H
#define INTERVAL_H

#include <iostream>
#include <string>
#include "base/FileParser.h"

class Interval
{
public:
Interval() {}

Interval(string &name, int startpt, int endpt) 
  : m_name(name), m_begin(startpt), m_end(endpt) {}

string Name() const {return m_name;}
int Begin() const {return m_begin;}
int End() const {return m_end;}

string GetName() const {return m_name;}
int GetStart() const {return m_begin;}
int GetEnd() const {return m_end;}

void SetName(string &name) {m_name=name;}
void SetBegin(int b) {m_begin=b;}
void SetEnd(int b) {m_end=b;}


int Size() const {return m_end-m_begin;}
bool Overlaps(Interval &e) 
{
  if (m_name != e.Name())
    return false;

  int beg = std::max(m_begin,e.Begin());
  int end = std::min(m_end,e.End());

  if (end-beg<0)
    return false;

  return true;
}

int AmtOverlap(Interval &e)
{
  if (m_name != e.Name())
    return 0;

  int beg = std::max(m_begin,e.Begin());
  int end = std::min(m_end,e.End());

  return (std::max(0,(end-beg)));

}


bool friend operator < (const Interval &lhs, const Interval &rhs) 
{
  if (lhs.m_name<rhs.m_name)
    return true;
  if (lhs.m_name>rhs.m_name)
    return false;

  if (lhs.m_begin<rhs.m_begin)
    return true;
  if (lhs.m_begin>rhs.m_begin)
    return false;

  if (lhs.m_end<rhs.m_end)
    return true;
  
  return false;

}

friend ostream & operator<< (ostream &ostrm, const Interval &e)
{
  ostrm << e.m_name <<"\t"<<e.m_begin<<"\t"<<e.m_end;
  return ostrm;
}


private:
string m_name;
int m_begin, m_end;


};




void LoadIntervals(string &filename, vector<Interval> &elements)
{
  elements.clear();
  FlatFileParser fp;
  fp.Open(filename);
  while (fp.ParseLine())
  {
    string name=fp.AsString(0);
    int start=fp.AsInt(1);
    int end = fp.AsInt(2);

    elements.push_back(Interval(name,start,end));
  }
}

#endif
