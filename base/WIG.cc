#include "base/WIG.h"

// WIG format
// track type... (ignore)
// variableStep chrom=chrN [span=windowSize] OR
// fixedStep chrom=chr1 start=1 step=200 [span=200]
// chromStartA dataValueA





void
ParseWig(const std::string &wigname, WIGFile* outfile)
{
  FlatFileParser f;
  StringParser sp;
  std::string line("");
  bool fixed_step = false;
  std::vector<char> seps;
  seps.push_back('=');
  float num_lines;
  FixedStepState fsstate;
  VariableStepState vsstate;
  std::vector<WIGLine> outvec;
  
  f.Open(wigname);
  num_lines = static_cast<float>(f.NumLines());
  DEBUGLOG(Stringify(num_lines) + " lines in WIG file.");  
  while (f.ParseLine()) {
    sp.SetLine(f.Line());
    if (Contains(sp.AsString(0), "fixedStep")) {      
      fixed_step = true;
      fsstate.chr = ChrToNum(sp.AsString(1).substr(6));
      fsstate.start = StringToInt(sp.AsString(2).substr(6));
      fsstate.step = StringToInt(sp.AsString(3).substr(5));
      if (sp.GetItemCount() > 4) {
        fsstate.span = StringToInt(sp.AsString(4).substr(5));
      }
      DEBUGLOG("Parsing chr " + Stringify(fsstate.chr));
    } else if (Contains(sp.AsString(0), "variableStep")) {
      fixed_step = false;
      vsstate.chr = ChrToNum(sp.AsString(1).substr(6));
      if (sp.GetItemCount() > 1) {
        vsstate.span = StringToInt(sp.AsString(1).substr(5));
      }
      DEBUGLOG("Parsing chr " + Stringify(vsstate.chr));
    }
    else {
      if (fixed_step) {
        ParseFixedStep(&sp, &fsstate, &outvec);
        outfile->lines[fsstate.chr].assign(outvec.begin(), outvec.end());
        outfile->span[fsstate.chr] = fsstate.span;
        outfile->step[fsstate.chr] = fsstate.step;
        outfile->start[fsstate.chr] = fsstate.start;
        outfile->is_fixed[fsstate.chr] = true;
      }
      else {
        ParseVariableStep(&sp, &vsstate, &outvec);
        outfile->lines[fsstate.chr].assign(outvec.begin(), outvec.end());
        outfile->span[fsstate.chr] = vsstate.span;
        outfile->is_fixed[fsstate.chr] = false;
      }
    }
  }
}

void 
ParseVariableStep(StringParser* sp, VariableStepState* state, std::vector<WIGLine>* outvec)
{
  WIGLine w;
  w.pos = sp->AsInt(0);
  w.val = sp->AsFloat(1);
  outvec->push_back(w);
}

void 
ParseFixedStep(StringParser* sp, FixedStepState* state, std::vector<WIGLine>* outvec)
{
  WIGLine w;
  w.pos = state->start + (state->count * state->step);
  w.val = sp->AsFloat(0);
  outvec->push_back(w);
  state->count++;  
}