#include "base/FileUtil.h"

bool FileExists(const std::string& fname)
{
  struct stat statbuffer;
  return (stat(fname.c_str(), &statbuffer) == 0);
}

inline void
GetAllFastaFilesInDirectory(const std::string& directory,
                            std::vector<std::string>* file_names)
{
}

inline std::string
GetFilenameWithoutPath(const std::string& filename)
{
    return std::string("");
}

inline std::string
GetFilenameWithoutExtension(const std::string& filename)
{
    return std::string("");
}



