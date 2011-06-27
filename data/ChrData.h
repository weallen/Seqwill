#ifndef CHRDATA_H_
#define CHRDATA_H_

#include <vector>
#include <string>
#include <ostream>

struct ChrData
{
  public:
  std::vector<std::string> chrnames;
  std::vector<int> chrlens;
  std::vector<std::string> tracknames;

  ChrData() {}
  virtual ~ChrData() {}
};

inline std::ostream& operator<< (std::ostream& s, const ChrData& c)
{
  s << "Chrs:" << std::endl;
  for (int i=0; i < c.chrnames.size(); ++i) {
    s << "Name: " << c.chrnames[i] << std::endl;
    s << "Len: " << c.chrlens[i] << std::endl;
  }
  s << "Tracks: " << std::endl;
  for (int i=0; i < c.tracknames.size(); ++i) {
    s << "Track: " << c.tracknames[i] << std::endl;
  }
}
#endif
