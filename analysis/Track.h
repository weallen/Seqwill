#ifndef TRACK_H_
#define TRACK_H_

#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace boost::numeric;

class Track
{
  public:
    Track(const string& name)
    : m_name(name) 
    {}

    virtual ~Track() {}
 
    const string& GetName() { return m_name; }

    ublas::vector<float> GetData() { return m_data; }
    void SetData(ublas::vector<float> data) { m_data = data; }

  private:
    string m_name;
    ublas::vector<float> m_data; 
};

#endif
