#ifndef WIG_H_
#define WIG_H_
#include <utility>
#include <vector>

#include "base/FileParser.h"
#include "base/StringUtil.h"

typedef std::map<std::string, std::string> KeyValMap;

struct WIGLine {
  int chr;
  int pos;
  float val;
};

void ParseWig(const std::string& wigname, std::vector<WIGLine>& out);

std::pair<std::string, std::string> ParseKeyValue(std::string kv);

#endif
