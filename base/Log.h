#ifndef LOG_H_
#define LOG_H_

#include <iostream>

const int INFO = 0;
const int WARNING = 1;
const int ERROR = 2;
const int FATAL = 3;

#define ERRORLOG(__Msg)                                                     \
    std::cerr << "Error: " << "[" << __func__ << ":" << __LINE__ << "] "    \
    << __Msg << std::endl

#define WARNLOG(__MSG)                                                      \
    std::cerr << "Warning: " << "[" << __func__ << ":" << __LINE__ << "] "  \
    << __Msg << std::endl

#endif
