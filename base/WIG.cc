#include "base/WIG.h"

// WIG format
// track type... (ignore)
// variableStep chrom=chrN [span=windowSize] OR
// fixedStep chrom=chr1 start=1 step=200 [span=200]
// chromStartA dataValueA

void
ParseWig(const std::string &wigname, std::vector<WIGLine> &out)
{
  FlatFileParser f;
  StringParser sp;
  std::string line("");
  bool fixed_step = false;
  std::vector<char> seps;
  seps.push_back('=');
  int start;
  int step;
  int span;
  std::string chrom;

  f.Open(wigname);
  while (f.ParseLine()) {
    sp.SetLine(f.Line());
    if (!sp.IsInt(0)) {
        if (Contains(sp.AsString(0), "fixedStep")) {
            fixed_step = true;
            chrom = ParseKeyValue(sp.AsString(1)).second;
            start = StringToInt(ParseKeyValue(sp.AsString(2)).second);
            step = StringToInt(ParseKeyValue(sp.AsString(3)).second);
            if (sp.GetItemCount() > 4) {
                span = StringToInt(ParseKeyValue(sp.AsString(4)).second);
            }
        } else if (Contains(sp.AsString(0), "variableStep")) {
          fixed_step = false;
          chrom = ParseKeyValue(sp.AsString(1)).second;
        if (sp.GetItemCount() > 1) {
          span = StringToInt(ParseKeyValue(sp.AsString(1)).second);
        }
      }
    }
  }
}

std::pair<std::string, std::string>
ParseKeyValue(std::string kv)
{
  static std::vector<char> seps;
  std::vector<std::string> toks;

  seps.clear();
  seps.push_back('=');
  Tokenize(kv, seps, toks);
  return std::make_pair(toks[0], toks[1]);
}
