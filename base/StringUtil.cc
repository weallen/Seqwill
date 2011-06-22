
#include "base/StringUtil.h"

int ExactPatternMatch(string orig, string pattern) {
    string::size_type pos = orig.find(pattern);
    if (pos == orig.npos) {
        return 0;
    }
    else {
        return 1;
    }
}

int IsWhitespace(char c) {
    return (c == ' ' or c == '\t' or c == '\n' or c == '\r');
}

int IsSpace(char c) {
    return (c == ' ' or c == '\t');
}

int ToWords(string &orig, vector<string> &words) {
    int curWordStart, curWordEnd;
    curWordStart = 0;
    while(curWordStart < orig.size()) {
        while (curWordStart < orig.size() and IsSpace(orig[curWordStart])) { 
            curWordStart++; 
        }
        curWordEnd = curWordStart;
        while (curWordEnd < orig.size() and !IsSpace(orig[curWordEnd])) { 
            curWordEnd++;
        }
        string word;
        if (curWordEnd != curWordStart) {
            word.assign(orig, curWordStart, curWordEnd - curWordStart);                     words.push_back(word);
        }
        curWordStart = curWordEnd;
    }
}


int AssignUntilFirstSpace(char *orig, int origLength, string &result) {
    int i;
    for (i = 0; i < origLength; i++ ){
        if (orig[i] == ' ' or orig[i] == '\t' or orig[i] == '\n' or orig[i] == '\r') {
            break;
        }
    }
    result.assign(orig, i);
    return i;
}


int PositionAfter(string& in, string& to_find, int startSearchAt)
{
  int max_pos = in.size() - to_find.size();
  for ( int i = startSearchAt; i <= max_pos; i++ )
  {
    int j;
    for ( j = 0; j < (int) to_find.size( ); ++j )
      if ( in[i + j] != to_find[j] )
        break;
    if ( j == (int) to_find.size( ) )
      return i;
  }
  return -1;
}


string After(string &s, string &t)
{
  string ret("");
  string::size_type pos = s.find(t);
  if ( pos != string::npos )
    ret = s.substr(pos+t.size(),s.size());

  return ret;
}


string Before(string &s, string &t)
{
  string ret("");
  string::size_type pos = s.find(t);
  if ( pos != string::npos )
    ret = s.substr(0,pos);

  return ret;

}


bool Contains(string &s, string &t)
{
  return (s.find(t) != string::npos);
}



string After(string &s, char *t)
{
  string ret("");
  string targ(t);
  string::size_type pos = s.find(targ);
  if ( pos != string::npos )
    ret = s.substr(pos+targ.size(),s.size());

  return ret;
}


string Before(string &s, char *t)
{
  string ret("");
  string::size_type pos = s.find(t);
  if ( pos != string::npos )
    ret = s.substr(0,pos);

  return ret;

}


bool Contains(string &s, char * t)
{

  return (s.find(t) != string::npos);

}


bool ContainsAt(string &s, string &t, int at)
{
  string::size_type pos = s.find(t);
  return ( (int) pos == at );
}



int Tokenize( const string &a_string,
	      vector<char> &separators,
	      vector<string> &tokens )
{

  tokens.clear();

  // Parse string.
  static vector< pair<int, int> > token_interval;
  token_interval.clear( );
  bool in_token = false;
  int token_start = 0;
  for (int ii=0; ii<(int)a_string.size( ); ii++) {
    bool matches_separator = false;
    for (int jj=0; jj<(int)separators.size( ); jj++)
      if ( separators[jj] ==a_string[ii] ) {
        matches_separator = true;
        break;
      }

    if ( matches_separator ) {
      if ( in_token ) {
        token_interval.push_back( pair<int, int>( token_start, ii ) );
        in_token = false;
      }
    }
    else {
      if ( !in_token ) {
        token_start = ii;
        in_token = true;
      }
    }
  }
  if ( in_token )    token_interval.push_back( pair<int, int>( token_start, (int)a_string.size() ) );

  // Create tokens.
  tokens.reserve( token_interval.size( ) );
  for (int ii=0; ii<(int)token_interval.size( ); ii++) {
    int token_length = token_interval[ii].second - token_interval[ii].first;
    tokens.push_back( a_string.substr( token_interval[ii].first, token_length ) );
  }

  // Return number of tokens.
  return (int)token_interval.size( );
}



int Tokenize( const string &a_string,
	      vector<string> &tokens )
{

  static vector<char> separators;
  separators.clear( );
  separators.push_back( ' ' );
  separators.push_back( '\t' );

  return Tokenize( a_string, separators, tokens );
}

