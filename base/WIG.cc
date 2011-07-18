#include "base/WIG.h"

// WIG format
// track type... (ignore)
// variableStep chrom=chrN [span=windowSize] OR
// fixedStep chrom=chr1 start=1 step=200 [span=200]
// chromStartA dataValueA



void
ParseWig(const std::string &wigname, std::vector<WIGLine>* outvec)
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
      std::cerr << "Chr: " << fsstate.chr << " , " << fsstate.start << " , " << fsstate.step << " , " << fsstate.span << std::endl;
    } else if (Contains(sp.AsString(0), "variableStep")) {
      fixed_step = false;
      vsstate.chr = ChrToNum(sp.AsString(1).substr(6));
      if (sp.GetItemCount() > 1) {
        vsstate.span = StringToInt(sp.AsString(1).substr(5));
      }
      DEBUGLOG("Parsing chr " + Stringify(vsstate.chr));
      std::cerr <<  "Chr: " << vsstate.chr << " , " << vsstate.span << std::endl;
    }
    else {
      if (fixed_step)
        ParseFixedStep(&sp, &fsstate, outvec);
      else
        ParseVariableStep(&sp, &vsstate, outvec);
    }
  }
}
void 
ParseVariableStep(StringParser* sp, VariableStepState* state, std::vector<WIGLine>* outvec)
{
  WIGLine w;
  w.pos = sp->AsInt(0);
  w.val = sp->AsFloat(1);
  w.chr = state->chr;
  outvec->push_back(w);
}

void 
ParseFixedStep(StringParser* sp, FixedStepState* state, std::vector<WIGLine>* outvec)
{
  WIGLine w;
  w.pos = state->start + (state->count * state->step);
  w.chr = state->chr;
  w.val = sp->AsFloat(0);
  outvec->push_back(w);
  state->count++;  
}

ChromosomeEnum 
ChrToNum(const std::string& chr) 
{
  if (chr.compare("chr1") == 0 || chr.compare("Chr1") == 0) 
    return kChr1;
  if (chr.compare("chr2")  == 0 || chr.compare("Chr2") == 0) 
    return kChr2;
  if (chr.compare("chr3")  == 0|| chr.compare("Chr3") == 0) 
    return kChr3;
  if (chr.compare("chr4") == 0 || chr.compare("Chr4")== 0) 
    return kChr4;
  if (chr.compare("chr5") == 0 || chr.compare("Chr5") == 0) 
    return kChr5;
  if (chr.compare("chr6") == 0 || chr.compare("Chr6") == 0) 
    return kChr6;
  if (chr.compare("chr7") == 0 || chr.compare("Chr7") == 0) 
    return kChr7;
  if (chr.compare("chr8") == 0 || chr.compare("Chr8") == 0) 
    return kChr8;
  if (chr.compare("chr9") == 0 || chr.compare("Chr9") == 0) 
    return kChr9;
  if (chr.compare("chr10") == 0 || chr.compare("Chr10") == 0) 
    return kChr10;
  if (chr.compare("chr11") == 0 || chr.compare("Chr11") == 0) 
    return kChr11;
  if (chr.compare("chr12") == 0 || chr.compare("Chr12") == 0) 
    return kChr12;
  if (chr.compare("chr13") == 0 || chr.compare("Chr13") == 0) 
    return kChr13;
  if (chr.compare("chr14") == 0 || chr.compare("Chr14") == 0) 
    return kChr14;
  if (chr.compare("chr15") == 0 || chr.compare("Chr15") == 0) 
    return kChr15;  
  if (chr.compare("chr16") == 0 || chr.compare("Chr16") == 0) 
    return kChr16;
  if (chr.compare("chr17") == 0 || chr.compare("Chr17") == 0) 
    return kChr17;
  if (chr.compare("chr18") == 0 || chr.compare("Chr18") == 0) 
    return kChr18;
  if (chr.compare("chr19") == 0 || chr.compare("Chr19") == 0) 
    return kChr19;
  if (chr.compare("chr20") == 0 || chr.compare("Chr20") == 0) 
    return kChr20;
  if (chr.compare("chrX") == 0 || chr.compare("ChrX") == 0) 
    return kChrX;
  if (chr.compare("chrY") == 0 || chr.compare("ChrY") == 0) 
    return kChrY;
  if (chr.compare("chrM") == 0 || chr.compare("ChrM") == 0) 
    return kChrM;
  return kChrUnknown;
} 
