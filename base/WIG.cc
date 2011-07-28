#include "base/WIG.h"

// WIG format
// track type... (ignore)
// variableStep chrom=chrN [span=windowSize] OR
// fixedStep chrom=chr1 start=1 step=200 [span=200]
// chromStartA dataValueA



void 
WIGParser::Open(const std::string& wigname)
{ 
  curr_state_ = new WIGState();
  filename_ = wigname; 
  fp_.Open(wigname); 
  num_lines_ = fp_.NumLines();
  DEBUGLOG(Stringify(fp_.NumLines()) + " lines in WIG file.");  
}


WIGLine
WIGParser::NextLine()
{
  curr_line_++;
  std::string line("");
  std::vector<char> seps;
  seps.push_back('=');
  bool chr_line = false;
  WIGLine w;
  fp_.ParseLine();
  sp_.SetLine(fp_.Line());
  if (Contains(sp_.AsString(0), "track")) {
    fp_.ParseLine();
    sp_.SetLine(fp_.Line()); // Read to next line
  }
  if (Contains(sp_.AsString(0), "fixedStep")) { 
    fixed_step_ = true;
    state_changed_ = true;
    chr_line = true;
    fsstate_.chr = ChrToNum(sp_.AsString(1).substr(6));
    fsstate_.start = StringToInt(sp_.AsString(2).substr(6));
    fsstate_.step = StringToInt(sp_.AsString(3).substr(5));
    if (sp_.GetItemCount() > 4) {
      fsstate_.span = StringToInt(sp_.AsString(4).substr(5));
    }
    fsstate_.count = 0;
    DEBUGLOG("Parsing chr " + Stringify(fsstate_.chr) + " at " + Stringify(fsstate_.step) + " from " + Stringify(fsstate_.start));
    delete curr_state_;
    curr_state_= new WIGState(fsstate_);
    
  } else if (Contains(sp_.AsString(0), "variableStep")) {
    fixed_step_ = false;
    chr_line = true;
    fixed_step_ = false;
    vsstate_.chr = ChrToNum(sp_.AsString(1).substr(6));
    if (sp_.GetItemCount() > 1) {
      vsstate_.span = StringToInt(sp_.AsString(1).substr(5));
    }
    DEBUGLOG("Parsing chr " + Stringify(vsstate_.chr));
    delete curr_state_;
    curr_state_ = new WIGState(vsstate_);
  }
  
  if (chr_line) {
    fp_.ParseLine();
    sp_.SetLine(fp_.Line()); // Read to next line
    state_changed_ = true;
  }
  if (!chr_line) {
    state_changed_ = false;
  }
  
  if (fixed_step_) {
    w.pos = fsstate_.start + (fsstate_.count * fsstate_.step);
    w.val = sp_.AsFloat(0);  
    fsstate_.count++;  
    return w;
  }
  else {
    w.pos = sp_.AsInt(0);
    w.val = sp_.AsFloat(1);
    return w;
  }
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
  if (chr.compare("chrX") == 0 || chr.compare("ChrX") == 0) 
    return kChrX;
  if (chr.compare("chrY") == 0 || chr.compare("ChrY") == 0) 
    return kChrY;
  if (chr.compare("chrM") == 0 || chr.compare("ChrM") == 0) 
    return kChrM;
  return kChrUnknown;
} 

 std::string ChrToString(ChromosomeEnum chr) 
 {
 if (chr == kChr1)
 return "chr1";
 if (chr == kChr2)
 return "chr2";
 if (chr == kChr3)
 return "chr3";
 if (chr == kChr4)
 return "chr4";
 if (chr == kChr5)
 return "chr5";
 if (chr == kChr6)
 return "chr6";
 if (chr == kChr7)
 return "chr7";
 if (chr == kChr8)
 return "chr8";
 if (chr == kChr9)
 return "chr9";
 if (chr == kChr10)
 return "chr10";
 if (chr == kChr11)
 return "chr11";
 if (chr == kChr12)
 return "chr12";
 if (chr == kChr13)
 return "chr13";
 if (chr == kChr14)
 return "chr14";
 if (chr == kChr15)
 return "chr15";
 if (chr == kChr16)
 return "chr16";
 if (chr == kChr17)
 return "chr17";
 if (chr == kChr18)
 return "chr18";
if (chr == kChr19)
  return "chr19";
 if (chr == kChrX)
 return "chrX";
 if (chr == kChrY)
 return "chrY";
 if (chr == kChrM)
 return "chrM";
 else
 return "unknown";
 }
 
