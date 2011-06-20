

#ifndef EST_
#define EST_


#include "aligns/AlignManager.h"
#include "base/Interval.h"

#include <map>


class EST
{
public:

EST() : mName("") 
{}

EST(string &n) 
  : mName(n), mDirection(""), mPartner("") {}

void SetDirection(string &d) { mDirection = d; }
void SetPartner(string &p) { mPartner = p; }
void AddInterval(interval &i) { mOnAssembly.push_back(i); }

void SetId(int id) { mId = id; }

string Direction() const  { return mDirection; }
string Partner() const { return mPartner; }
string Name() const { return mName; }
int Id() const { return mId; }

bool Is5primeEnd() { return (mDirection == "F"); }
bool Is3primeEnd() { return (mDirection == "R"); } 

bool IsEmpty() { return (mName.empty()); }

bool SpansSupers() const { return (mOnAssembly.size() > 1); }
vector<interval> & Intervals()  { return mOnAssembly; }
int NumPlacements() const { return mOnAssembly.size(); }

const interval & GetInterval(int indx) const
{ 
  if  ( indx < 0 || indx >= (int) mOnAssembly.size() )
  {
    cout << "invalid index: " << indx <<"; number of intervals=" << mOnAssembly.size() << endl;
    exit(-1);
  }
  return mOnAssembly[indx]; 
}


bool Overlaps(const EST &e);


private:

string mName, mDirection, mPartner;
vector<interval> mOnAssembly;
int mId;


};


bool EST::Overlaps(const EST &e)
{
  for ( int i=0; i<(int) mOnAssembly.size(); ++i )
  {
    string superId = mOnAssembly[i].ReferenceName(); 

    int numPlacements = e.NumPlacements();
    for ( int j=0; j<numPlacements; ++j )
    {
      if ( e.GetInterval(j).ReferenceName() != superId )
	continue;

      if ( mOnAssembly[i].HasOverlap(e.GetInterval(j) ) )
	return true;
    }
  }

  return false;
    
}




class EST_Mgr : public SpinesAlignsMgr
{

public:
EST_Mgr(string &dir) : SpinesAlignsMgr(dir) {  }


typedef pair<EST,EST> cDNA;


void SetupMap();
void CollectESTs();

vector<EST> GetESTs() { return m_ESTs; }
EST  GetPartnerEST(const EST &e);

bool OverlapsPartner(const EST &e);

void Form_cDNAs();
vector<cDNA> GetcDNAs() { return m_cDNAs; }

int NumcDNAs() { return m_cDNAs.size(); }


// only one end of the cDNA either a) exists or b) aligns
int Num5pOnly();
int Num3pOnly();
int NumBothEnds();


int Overlapping_Partner_Count();

int Overlapping_cDNA_Count();
int Overlapping_5prime_Count();
int Overlapping_3prime_Count();


int NumClusters() { return mClusters.size(); }
vector< vector<cDNA> > Clusters() { return mClusters; }

private:

typedef multimap<string,QCMark *> estMap;
typedef estMap::iterator estMapIter;
estMap m_estMap;

map<string,int> mNameToIndex;

vector<EST> m_ESTs;
vector<cDNA> m_cDNAs;
vector< vector<cDNA> > mClusters;

bool Overlapping_cDNAs(cDNA &c1,cDNA &c2);

};


int EST_Mgr::Num5pOnly()
{
  int count(0);
  if ( m_cDNAs.empty() )
    return count;

  for ( int i=0; i<(int) m_cDNAs.size(); ++i )
  {
    cDNA & cdna = m_cDNAs[i];

    bool have5p(false), have3p(false);
    if ( cdna.first.Is5primeEnd() && !cdna.first.IsEmpty() )
      have5p=true;
    if ( !have5p && (cdna.second.Is5primeEnd() && !cdna.second.IsEmpty()) )
      have5p=true;

    if ( cdna.first.Is3primeEnd() && !cdna.first.IsEmpty() )
      have3p=true;
    if ( !have3p && (cdna.second.Is3primeEnd() && !cdna.second.IsEmpty()) )
      have3p=true;
    
    if ( have5p && !have3p )
      ++count;
  }
      
  return count;

}

int EST_Mgr::Num3pOnly()
{
  int count(0);
  if ( m_cDNAs.empty() )
    return count;

  for ( int i=0; i<(int) m_cDNAs.size(); ++i )
  {
    cDNA & cdna = m_cDNAs[i];

    bool have5p(false), have3p(false);
    if ( cdna.first.Is5primeEnd() && !cdna.first.IsEmpty() )
      have5p=true;
    if ( !have5p && (cdna.second.Is5primeEnd() && !cdna.second.IsEmpty()) )
      have5p=true;

    if ( cdna.first.Is3primeEnd() && !cdna.first.IsEmpty() )
      have3p=true;
    if ( !have3p && (cdna.second.Is3primeEnd() && !cdna.second.IsEmpty()) )
      have3p=true;
    
    if ( !have5p && have3p )
      ++count;
  }
      
  return count;

}

int EST_Mgr::NumBothEnds()
{
  int count(0);
  if ( m_cDNAs.empty() )
    return count;

  for ( int i=0; i<(int) m_cDNAs.size(); ++i )
  {
    cDNA & cdna = m_cDNAs[i];

    bool have5p(false), have3p(false);
    if ( cdna.first.Is5primeEnd() && !cdna.first.IsEmpty() )
      have5p=true;
    if ( !have5p && (cdna.second.Is5primeEnd() && !cdna.second.IsEmpty()) )
      have5p=true;

    if ( cdna.first.Is3primeEnd() && !cdna.first.IsEmpty() )
      have3p=true;
    if ( !have3p && (cdna.second.Is3primeEnd() && !cdna.second.IsEmpty()) )
      have3p=true;
    
    if ( have5p && have3p )
    {
      ++count;
      // 
      if ( !OverlapsPartner(cdna.first) )
      {
	if  ( cdna.first.SpansSupers() ||
	      cdna.second.SpansSupers() )
	  continue;
	
	interval v1 = cdna.first.GetInterval(0);
	interval v2 = cdna.second.GetInterval(0);

	if ( v1.ReferenceName() != v2.ReferenceName() )
	  continue;

	cout << cdna.first.Name() 
	     << endl;
	cout << cdna.second.Name()
	     << endl;
      }
    }
  }
      
  return count;

}


int EST_Mgr::Overlapping_Partner_Count()
{
  int count(0);
  if ( m_cDNAs.empty() )
    return count;

  for ( int i=0; i<(int) m_cDNAs.size(); ++i )
  {
    
    if ( OverlapsPartner(m_cDNAs[i].first) )
      ++count;
  }
      
  return count;

}


int EST_Mgr::Overlapping_5prime_Count()
{
  if ( m_cDNAs.empty() )
    return 0;

  vector<int> used(m_cDNAs.size(),0);
  int count(0);

  for ( int i=0; i<(int) m_cDNAs.size(); ++i)
  {
    if ( used[i] )
      continue;

    cDNA &c = m_cDNAs[i];

    EST this_est = c.first;
    if ( c.first.Is3primeEnd() )
      this_est = c.second;

    if ( this_est.IsEmpty() ) 
    {
      used[i]=1;
      continue;
    }

    for ( int j=i+1; j<(int) m_cDNAs.size(); ++j )
    {
      if ( used[j] )
	continue;
      
      EST other_est = m_cDNAs[j].first;
      if ( other_est.Is3primeEnd() )
	other_est = m_cDNAs[j].second;

      if ( other_est.IsEmpty() )
      {
	used[j]=1;
	continue;
      }

      if ( this_est.Overlaps(other_est) )
      {
	++count;
	used[j]=1;
      }

    }

    used[i]=1;
  }

  return count;
}

int EST_Mgr::Overlapping_3prime_Count()
{
  if ( m_cDNAs.empty() )
    return 0;

  vector<int> used(m_cDNAs.size(),0);
  int count(0);

  for ( int i=0; i<(int) m_cDNAs.size(); ++i)
  {
    if ( used[i] )
      continue;

    cDNA &c = m_cDNAs[i];

    EST this_est = c.first;
    if ( c.first.Is5primeEnd() )
      this_est = c.second;

    if ( this_est.IsEmpty() ) 
    {
      used[i]=1;
      continue;
    }

    for ( int j=i+1; j<(int) m_cDNAs.size(); ++j )
    {
      if ( used[j] )
	continue;
      
      EST other_est = m_cDNAs[j].first;
      if ( other_est.Is5primeEnd() )
	other_est = m_cDNAs[j].second;

      if ( other_est.IsEmpty() )
      {
	used[j]=1;
	continue;
      }

      if ( this_est.Overlaps(other_est) )
      {
	++count;
	used[j]=1;
      }

    }

    used[i]=1;
  }

  return count;
}



int EST_Mgr::Overlapping_cDNA_Count()
{
  int count(0);
  vector<int> used(m_cDNAs.size(),0);

  for ( int i=0; i<(int) m_cDNAs.size(); ++i)
  {
    if ( used[i] )
      continue;

    vector<cDNA> this_cluster;
    this_cluster.push_back(m_cDNAs[i]);
   
    for ( int j=i+1; j<(int) m_cDNAs.size(); ++j )
    {
      if ( used[j] )
	continue;
      
      if ( Overlapping_cDNAs(m_cDNAs[i],m_cDNAs[j]) )
      {
	++count;
	used[j]=1;
	this_cluster.push_back(m_cDNAs[j]);
      }
    }
    
    used[i]=1;
    
    if ( this_cluster.size() > 1 )
    {
      ++count;
      mClusters.push_back(this_cluster);
    }
  }
  return count;
}


bool  EST_Mgr::Overlapping_cDNAs(cDNA &c1,cDNA &c2)
{
  return ( c1.first.Overlaps(c2.first) ||
	   c1.second.Overlaps(c2.first) ||
	   c1.first.Overlaps(c2.second ) ||
	   c1.second.Overlaps(c2.second) );
}


void EST_Mgr::Form_cDNAs()
{
  vector<int> used(m_ESTs.size(),0);
  EST empty;
  
  for ( int i=0;i<(int) m_ESTs.size(); ++i )
  {
    if ( used[i] )
      continue;

    EST &e = m_ESTs[i];

    string partner = m_ESTs[i].Partner();
    map<string,int>::iterator mIter = mNameToIndex.find(partner);
    if ( mIter == mNameToIndex.end() )
    {
      m_cDNAs.push_back(make_pair(e,empty));
    }
    else
    {
      EST &p = m_ESTs[mIter->second];
      m_cDNAs.push_back(make_pair(e,p));
      used[mIter->second]=1;
    }
    
    used[i]=1;
  }
}


bool EST_Mgr::OverlapsPartner(const EST &e)
{
 EST p = GetPartnerEST(e);
  if ( p.IsEmpty() )
    return false;
  
  
  return ( e.Overlaps(p) );

}

EST EST_Mgr::GetPartnerEST(const EST &e)
{
  string partner = e.Partner();
  map<string,int>::iterator mIter = mNameToIndex.find(partner);
  if ( mIter == mNameToIndex.end() )
  {
    return (EST());
  }

  return (m_ESTs[mIter->second]);

}



void EST_Mgr::SetupMap()
{
  int num_targets = m_details.GetNumSupers();

//   cout << "in setup: " << num_targets << endl;
  for ( int t=0; t<num_targets; ++t )
  {
    int num_hits_this_target(m_details.GetCountForSuper(t));

    // no alignment to this target super
    if ( num_hits_this_target == 0 )
      continue;
    
    for ( int q=0; q<num_hits_this_target; ++q )
    {
      QCMark & mark = m_details.MarkForSuper(t,q);
      if ( mark.GetType() != "aligned" )
	continue;
      
      int id = atoi(QueryNameFromMark(mark).c_str());
      string name = mMetainfoMgr->GetNameFromId(id);
      m_estMap.insert(make_pair(name,&mark));
    }
  }
}

void EST_Mgr::CollectESTs()
{

  const ReadMetaInfoBroker * metainfo = mMetainfoMgr->FullAccess();

  set<int> alreadyListed;
 
  int numQ = m_details.GetNumSupers();
  for ( int i=0; i<numQ; ++i )
  {
    int count = m_details.GetCountForSuper(i);
    if ( count == 0 )
      continue;

    for ( int j=0; j<count; ++j )
    {

      QCMark & mark = m_details.MarkForSuper(i,j);
      if ( mark.GetType() != "aligned" )
	continue;
      
      int id = atoi(QueryNameFromMark(mark).c_str());      
      string est_name(mMetainfoMgr->GetNameFromId(id));

      if ( !alreadyListed.insert(id).second )
	continue;

      pair<estMapIter,estMapIter> range = m_estMap.equal_range(est_name);
      if ( range.first==range.second )
	continue;
    

      EST this_est(est_name);
      string direction(metainfo->Direction(est_name));
      this_est.SetDirection(direction);
      this_est.SetId(id);
      
      string partner(mMetainfoMgr->GetPartner(est_name));
      this_est.SetPartner(partner);
      
      int thisSuper=range.first->second->GetSuper();
      int start=range.first->second->GetStartOnSuper();
      int end=range.first->second->GetEndOnSuper();
      
      string ref(Stringify(thisSuper));
      string ivid("1");
      interval sI(ref,ivid,start,end);
      
      ++range.first;
      for (; range.first != range.second; ++range.first )
      {
	
	if ( range.first->second->GetSuper() == thisSuper )
	{
	  start=std::min(start,range.first->second->GetStartOnSuper());
	  end=std::max(end,range.first->second->GetEndOnSuper());
	  
	  sI.SetBegin(start);
	  sI.SetEnd(end);
	  
	}
	else
	{
	  this_est.AddInterval(sI);
	  // 	m_ESTs.push_back(this_est);
	  
	  thisSuper = range.first->second->GetSuper();
	  start = range.first->second->GetStartOnSuper();
	  end = range.first->second->GetEndOnSuper();
	  string sup(Stringify(thisSuper));
	  sI.SetReferenceName(sup);
	  sI.SetBegin(start);
	  sI.SetEnd(end);
	}
      }
      
      this_est.AddInterval(sI);
      m_ESTs.push_back(this_est);

    }

  }

  for ( int i=0; i<(int) m_ESTs.size(); ++i )
    mNameToIndex.insert(make_pair(m_ESTs[i].Name(),i));

}



#endif
