#ifndef CHRDATA_H_
#define CHRDATA_H_

#include <vector>
#include <string>

class ChrData
{
  public:
  std::vector<std::string> chrnames;
  std::vector<int> chrlens;
  std::vector<std::string> tracknames;

  ChrData() {}
  virtual ~ChrData() {}

};
#endif
