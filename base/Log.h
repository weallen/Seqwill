#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <stdio.h>
#include <stdarg.h>

const int INFO = 0;
const int WARNING = 1;
const int ERROR = 2;
const int FATAL = 3;

#define ERRORLOG(__Msg)                                                     \
    std::cerr << "Error: " << "[" << __FILE__  << ":" << __func__ << ":" << __LINE__ << "] "    \
    << __Msg << std::endl

#define WARNLOG(__Msg)                                                      \
    std::cerr << "Warning: " << "[" << __FILE__ << ":" << __func__ << ":" << __LINE__ << "] "  \
    << __Msg << std::endl

#define DEBUGLOG(__Msg)\
    std::cout << "Debug: " << "[" << __FILE__ << ":" << __func__ << ":" << __LINE__ << "] " \
    << __Msg << std::endl;\

#endif
