//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif



#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{
  commandArg<string> iStringCmmd("-i","input fasta file");
  commandArg<string> cStringCmmd("-c","string to search for, e.g. GATTACA","");
  commandArg<bool> rStringCmmd("-repeats","count repeat content",false);
  commandArg<bool> gStringCmmd("-gc","count GC content",false);
  commandArg<int> wStringCmmd("-gcByWindow","window size for GC content",-1);
  commandArg<int> rwStringCmmd("-repeatByWindow","window size for repeat content",-1);


  commandLineParser P(argc,argv);
  P.SetDescription("Counts occurrences of a particular string in a fasta file without overlap.");
  P.registerArg(iStringCmmd);
  P.registerArg(cStringCmmd);
  P.registerArg(rStringCmmd);
  P.registerArg(gStringCmmd);
  P.registerArg(rwStringCmmd);
  P.registerArg(wStringCmmd);

  P.parse();

  string fasta = P.GetStringValueFor(iStringCmmd);
  string pattern = P.GetStringValueFor(cStringCmmd);
  bool repeat = P.GetBoolValueFor(rStringCmmd);
  bool gc = P.GetBoolValueFor(gStringCmmd);
  int rw = P.GetIntValueFor(rwStringCmmd);
  int w = P.GetIntValueFor(wStringCmmd);

  while(pattern[0]=='N') pattern.erase(0,1);
  while(pattern[pattern.length()-1]=='N') pattern.erase(pattern.length()-1,1);

  FlatFileParser fastaparse;
  fastaparse.Open(fasta);
  
  string templine; 
  int state;
  int count=0;
  int pos=0;
  int lastend;
  string name;
  bool start;
 
  int x1=0;
  int x2=0;
  int x3=0;
  int x4=0;
  int x5=0;

  if(w>0) x1=1;
  if(pattern.size()>0) x2=1;
  if(repeat) x3=1;
  if(gc) x4=1;
  if(rw>0) x5=1;

  if(x1+x2+x3+x4+x5 != 1) {
    cerr << "Choose one: -c , -repeats , -gc , -gcbyWindow , -repeatByWindow" << endl;
    exit(-1);
  }


  if(w > 0) {

    cout << "SequenceName\tstart\tend\tA+Tcontent\tG+Ccontent\tGC%" << endl;
    
    start = 0;
    int gcGC = 0;
    int atAT = 0;
    lastend = 0;

    while(fastaparse.GetLine(templine)) {
      
      if(templine[0]=='>') {
	if(start) cout << name << "\t" << lastend << "\t" << pos-1 << "\t" <<  atAT << "\t" << gcGC << "\t" << (double)gcGC/((double)gcGC+(double)atAT) << endl;
	templine.erase(0,1);
	name=templine;
	gcGC = 0;
	atAT = 0;
	pos = 0;
	start = 0;
	lastend = 0;
	continue;
      }

      while(!templine.empty()) {

        if(pos % w == 0) {
          if(start) cout << name << "\t" << pos-w << "\t" << pos-1 << "\t" << atAT << "\t" << gcGC << "\t" << (double)gcGC/((double)gcGC+(double)atAT)  << endl;gcGC=0;
          atAT=0;
	  lastend=pos;
        }

	start = 1;
	if(toupper(templine[0])=='G' || toupper(templine[0])=='C') gcGC++;
	else if(toupper(templine[0])=='A' || toupper(templine[0])=='T') atAT++;
      
	templine.erase(0,1);
        pos++;

	

      }
    }
    
    cout << name << "\t" << lastend << "\t" << pos-1 << "\t" <<  atAT << "\t" << gcGC << "\t" << (double)gcGC/((double)gcGC+(double)atAT) << endl;
  
    return 0;

  }



if(rw > 0) {

    cout << "SequenceName\tstart\tend\tNonRepeats\tRepeats\tRepeat%" << endl;
    
    start = 0;
    int actgn = 0;
    int ACTGN = 0;
    lastend = 0;

    while(fastaparse.GetLine(templine)) {
      
      if(templine[0]=='>') {
	if(start) cout << name << "\t" << lastend << "\t" << pos-1 << "\t" <<  ACTGN << "\t" << actgn << "\t" << (double)actgn/((double)actgn+(double)ACTGN) << endl;
	templine.erase(0,1);
	name=templine;
	actgn = 0;
	ACTGN = 0;
	pos = 0;
	start = 0;
	lastend = 0;
	continue;
      }

      while(!templine.empty()) {

        if(pos % rw == 0) {
          if(start) cout << name << "\t" << pos-rw << "\t" << pos-1 << "\t" << ACTGN << "\t" << actgn << "\t" << (double)actgn/((double)actgn+(double)ACTGN)  << endl; 
	  actgn=0;
          ACTGN=0;
	  lastend=pos;
        }

	start = 1;
	if(islower(templine[0])) actgn++;
	else if(isupper(templine[0])) ACTGN++;
      
	templine.erase(0,1);
        pos++;

	

      }
    }
    
    cout << name << "\t" << lastend << "\t" << pos-1 << "\t" <<  ACTGN << "\t" << actgn << "\t" << (double)actgn/((double)actgn+(double)ACTGN) << endl;
  
    return 0;

  }




  if(!repeat && !gc) {

    state = 0;
    start = 0;
    count = 0;
    cout << "SequenceName\tCount_WithoutOverlaps" << endl;

    while(fastaparse.GetLine(templine))
      {

	if(templine[0]=='>') {
	  if(start) cout << name << "\t" << count << endl;  
	  templine.erase(0,1);
	  name=templine;
	  pos = 0;
	  start = 1;
	  count = 0;
	  state = 0;
	  continue;
	}



        while(!templine.empty())
          {
	   
	    if(toupper(templine[0])==toupper(pattern[state]) || toupper(pattern[state]=='N')) state++;
	      else state=0;

	      if(state==pattern.length())
	      {
                state=0;
	        count++;
	      }
	    
	  
	    templine.erase(0,1);
	    pos++;
	  }
      }
      
    cout << name << "\t" << count << endl;  
    return 0;
  } 



  if(repeat) {
    int actgn = 0;
    int ACTGN = 0;
    start = 0;

    cout << "SequenceName\tRepeats\tNon-Repeats\tRepeat%" << endl;

    while(fastaparse.GetLine(templine))
    {
    
	if(templine[0]=='>') {
	  if(start) cout << name << "\t" << actgn << "\t" << ACTGN << "\t" << (double)actgn/((double)actgn + (double)ACTGN) << endl;
	  templine.erase(0,1);
	  name=templine;
	  pos = 0;
	  actgn = 0;
	  ACTGN = 0;
	  start = 1;
	  continue;
	}



      while(!templine.empty())
        {
          

	    if(templine[0]=='a' || templine[0]=='c' || templine[0]=='t' || templine[0]=='g' || templine[0]=='n') actgn++;
            else if(templine[0]=='A' || templine[0]=='C' || templine[0]=='T' || templine[0]=='G' || templine[0]=='N') ACTGN++; 	  

	  
	  templine.erase(0,1);
          pos++;
        
	}
    }
      
    cout << name << "\t" << actgn << "\t" << ACTGN << "\t" << (double)actgn/((double)actgn + (double)ACTGN) << endl;

    return 0;
  }


  if(gc) {
    start = 0;
    int gcGC = 0;
    int atAT = 0;

    cout << "SequenceName\tA+Tcontent\tG+Ccontent\tGC%" << endl;

    while(fastaparse.GetLine(templine))
      {

	if(templine[0]=='>') {
	  if(start) cout << name << "\t" << atAT << "\t" << gcGC << "\t" << (double)gcGC/((double)gcGC+(double)atAT) << endl;
	  templine.erase(0,1);
	  name=templine;
	  pos = 0;
	  gcGC = 0;
	  atAT = 0;
	  start = 1;
	  continue;
	}



	while(!templine.empty())
	  {
          

	      if(templine[0]=='g' || templine[0]=='c' || templine[0]=='G' || templine[0]=='C') gcGC++;
	      else if(templine[0]=='a' || templine[0]=='t' || templine[0]=='A' || templine[0]=='T') atAT++; 	  

	  
	    templine.erase(0,1);
	    pos++;
        
	  }
      }
      
    cout << name << "\t" << atAT << "\t" << gcGC << "\t" << (double)gcGC/((double)gcGC+(double)atAT) << endl;

      return 0;
  }


  return 0;
}
    
