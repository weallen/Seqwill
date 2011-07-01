#include "io/TrackReader.h"

template<typename T>
void TrackReader<T>::Open(const char* dirname)
{
  std::string s(dirname);
  Open(s);
}

template<typename T>
void TrackReader<T>::Open(const std::string& dirname)
{
}

template<typename T>
void TrackReader<T>::Close()
{
  Flush();
}

// Deallocate everything
template<typename T>
void TrackReader<T>::Flush()
{
}

