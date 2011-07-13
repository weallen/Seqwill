// Copyright (c) 2010, Pacific Biosciences of California, Inc.
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright notice, thislist of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//     * Neither the name of Pacific Biosciences nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. I
//   N NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS BE LIABLEFRAY//   // DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDIN
//   G, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWA

#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdint.h>

using namespace std;

// Stuff from PacBio

int ExactPatternMatch(string orig, string pattern);
int IsWhitespace(char c);
int IsSpace(char c);
int ToWords(string& orig, vector<string>& words);
int AssignUntilFirstSpace(char* orig, int origLength, string& result);


// Stuff from Spines
string After(string &s, string &t);
string Before(string &s, string &t);
bool Contains(string &s, string &t);

string After(const string &s, const char *t);
string Before(const string &s, const char *t);
bool Contains(const string &s, const char *t);

bool ContainsAt(string &s, string &t, int at);

int PositionAfter(string &in, string& s, int startSearchAt);

inline string Stringify(int x)
{
  ostringstream out;
  out << x;
  return (out.str());
}

inline int StringToInt(string __s) 
{
  return atoi(__s.c_str());
}

int Tokenize( const string &a_string,
	      vector<char> &separators,
	      vector<string> &tokens );


int Tokenize( const string &a_string,
	      vector<string> &tokens );


  inline uint64_t HashString(const char* __s)
  {
    uint64_t hash = 0xcbf29ce484222325ull;
    for ( ; *__s; ++__s) {
      hash *= 1099511628211ull;
      hash ^= *__s;
    }
    return hash;
  }


#endif
