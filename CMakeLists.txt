option(build_tests "Build tests." ON)
  
cmake_minimum_required (VERSION 2.0)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_C_COMPILER "/usr/local/bin/clang")

set(CMAKE_CXX_COMPILER "/usr/local/bin/clang++")
project (mediphmm)

#set(CMAKE_CXX_FLAGS "-Wall -Wshadow")
#set(cxx_exception_flags "-fexceptions")
#set(cxx_no_exception_flags "-fnoexceptions")

set(EXECUTABLE_OUTPUT_PATH bin)

include_directories(".")
include_directories("/usr/local/include/bamtools")

# include Google Test
set(gtest_SOURCE_DIR "~/src/googletest-read-only")
include_directories("${gtest_SOURCE_DIR}/include" ${gtest_SOURCE_DIR})
link_directories("${gtest_SOURCE_DIR}")

# libraries
add_library(util util/SComm.cc util/SysTime.cc util/mutil.cc base/ErrorHandling.cc base/FileParser.cc base/StringUtil.cc) 

add_library(analysis analysis/Genome.cc analysis/Track.cc analysis/Chromosome.cc analysis/Bam.cc)
target_link_libraries(analysis hdf5_cpp hdf5 bamtools util)

# binaries
add_executable(LoadSeq tool/LoadSeq.cc)
add_executable(LoadData tool/LoadData.cc)
target_link_libraries(LoadSeq util analysis)
target_link_libraries(LoadData util analysis)

# testing
function(cxx_test name)
  add_executable(${name} "test/${name}.cc")
  target_link_libraries(${name} analysis util gtest)
endfunction()

cxx_test(GenomeTest)
cxx_test(ChromosomeTest)
