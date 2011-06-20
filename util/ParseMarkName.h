////////////
#ifndef _PARSEMARKNAME_H_
#define _PARSEMARKNAME_H_

#include "base/StringUtil.h"
#include <stdlib.h>
#include <string.h>

// Stupid and inefficient function to parse a mark name separated by ':'

inline void ParseMarkName(vector<string> & tokens, const string & name) 
{
  const char * p = name.c_str();
 
  int len = strlen(p);

  string all;
  for (int i=0; i<len; i++) {
    if (p[i] == ':' || p[i] == ';') {
      tokens.push_back(all);      
      all = "";
      continue;
    }
    all += p[i];
  }
  tokens.push_back(all);
} 

// Fill misjoin, indel, and errors with the proper floats (range [0,100)).

inline void ParseMarkNameFloat( const string &name,
				double &misjoin,
				double &indel,
				double &errors )
{
  vector<string> tokens;
  vector<char> sep;
  sep.push_back( ':' );
  sep.push_back( ';' );
  Tokenize( name, sep, tokens );

  //ForceAssert( tokens.size( ) == 6 );
  //ForceAssert( tokens[0] == "p_misjoin" );
  //ForceAssert( tokens[2] == "p_indel" );
  //ForceAssert( tokens[4] == "p_baseerrors" );

  misjoin = 100.0 * atof(tokens[1].c_str());
  indel = 100.0 * atof(tokens[3].c_str());
  errors = 100.0 * atof(tokens[5].c_str());
}

#endif //_PARSE_MARK_NAME_H_

