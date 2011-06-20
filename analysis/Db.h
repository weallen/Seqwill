#ifndef DB_H_
#define DB_H_

#include <iostream>
#include <string>
#include <sqlite3.h>

using namespace std;

class DbMgr
{
  public:
    Db(const string& dirname) 
      : m_dbname(dirname + "/data.db")
    { }
    virtual ~Db() {}

  
    int Open() {
      Close()
      return sqlite3_open(m_dbname.c_str(), &m_db);   
    }
    
    void Close() {
      if (m_db) {
        sqlite3_close(m_db);
        m_db = NULL;
      }
    }
  
    // E.g. Execute("CREATE TABLE IF NOT EXISTS blah ( a INTEGER )"
    int Execute(const string& query) {
      return sqlite3_exec(m_db, query.c_str(), NULL, NULL, NULL);
    }

    void 
  private:
    string m_dbname;
    sqlite3* m_db;
};

#endif

