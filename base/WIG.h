#ifndef WIG_H_
#define WIG_H_

#include <iostream>
#include <utility>
#include <vector>
#include <map>
#include <limits>
#include <stdio.h>

#include "base/FileParser.h"
#include "base/StringUtil.h"
#include "base/Types.h"
#include "base/Log.h"

struct WIGLine;

struct WIGFile {
  std::map<ChromosomeEnum, std::vector<WIGLine> > lines;
  std::map<ChromosomeEnum, int> span;
  std::map<ChromosomeEnum, int> step;
  std::map<ChromosomeEnum, int> start;
  std::map<ChromosomeEnum, bool> is_fixed;
};

struct WIGLine {
    int pos;
    float val;
};

struct FixedStepState {
  FixedStepState()
  :   chr(kChrUnknown) 
  , start(-1)
  , step(-1)
  , span(-1)
  , count(0)
  {}
  
  FixedStepState(ChromosomeEnum ichr,
                 int istart, int istep,
                 int ispan)
  { chr = ichr; start = istart; 
    step = istep; span = ispan; }
  
  ChromosomeEnum chr;
  int start;
  int step;
  int span;
  int count;
};


struct VariableStepState {
  VariableStepState()
  : chr(kChrUnknown)
  , span(-1)
  {}
  
  VariableStepState(ChromosomeEnum ichr, 
                    int ispan)
  { chr = ichr; span = ispan; }
  
  ChromosomeEnum chr;
  int span;
};

void ParseWig(const std::string& wigname, WIGFile* outfile);  

void ParseVariableStep( StringParser* f,  VariableStepState* state, std::vector<WIGLine>* outvec);

void ParseFixedStep( StringParser* f,  FixedStepState* state, std::vector<WIGLine>* outvec);

#endif
