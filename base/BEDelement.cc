#include "base/BEDelement.h"

void LoadBEDfile(string &filename, vector<BEDelement> &elements)
{
  elements.clear();
  FlatFileParser fp;

  fp.Open(filename);
  while (fp.ParseLine())
  {
    string name=fp.AsString(0);
    int start=fp.AsInt(1);
    int end = fp.AsInt(2);
    string genename = ""; //fp.AsString(3);

    elements.push_back(BEDelement(name,start,end,genename));
  }
}
