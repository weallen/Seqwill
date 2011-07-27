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
struct WIGState {
  
  WIGState() 
  : chr(kChrUnknown)
  , span(-1)
  , step(-1)
  , start(-1)
  , is_fixed(false)
  {}
  
  WIGState(ChromosomeEnum ichr, int ispan,
           int istep, int istart, bool iis_fixed)
  { chr = ichr; span = ispan; step = istep; start = istart;
    is_fixed = iis_fixed; }
  
  WIGState(const FixedStepState& s)
  { chr = s.chr; span = s.span; step = s.step;
    start = s.start; is_fixed = true; }
  
  WIGState(const VariableStepState& s)
  { chr = s.chr, span = s.span; is_fixed = false;
    start = -1; step = -1; }
  
  ChromosomeEnum chr;
  int span;
  int step;
  int start;
  bool is_fixed;
};

class WIGParser 
{
public:
  WIGParser() 
  : state_changed_(false)
  , num_lines_(0)  
  , filename_("")
  , curr_line_(0)
  , fixed_step_(false)
  {}
  
  virtual ~WIGParser() { }  
  void Open(const std::string& wigname);
  bool HasNextLine() { return curr_line_ < num_lines_;}
  WIGLine NextLine();  

  bool StateChanged() const { return state_changed_; }
  WIGState curr_state() { return *curr_state_; }
  int num_lines() const { return num_lines_; }
private:  
  WIGLine ParseVariableStep();
  WIGLine ParseFixedStep();
    
  bool state_changed_;
  int num_lines_;
  std::string filename_;
  FlatFileParser fp_;
  FixedStepState fsstate_;
  VariableStepState vsstate_;
  StringParser sp_;
  WIGState* curr_state_;
  int curr_line_;
  bool fixed_step_;
};

ChromosomeEnum 
ChrToNum(const std::string& chr);

std::string ChrToString(ChromosomeEnum chr);

#endif
