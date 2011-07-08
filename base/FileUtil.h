#ifndef FILEUTIL_H_
#define FILEUTIL_H_

#include <iostream>
#include <vector>
#include <sys/stat.h>

bool FileExists(const std::string& fname);

inline void
GetAllFastaFilesInDirectory(const std::string& directory,
                            std::vector<std::string>* file_names);

inline std::string
GetFilenameWithoutPath(const std::string& filename);

inline std::string
GetFilenameWithoutExtension(const std::string& filename);

#endif
