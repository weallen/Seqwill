cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_C_COMPILER "clang")

set(CMAKE_CXX_COMPILER "clang++")
project(Seqwill)

set(CMAKE_CXX_FLAGS "-g")
#set(cxx_exception_flags "-fexceptions")
#set(cxx_no_exception_flags "-fnoexceptions")

set(EXECUTABLE_OUTPUT_PATH bin)

include_directories(".")
include_directories("/usr/local/include/bamtools")
include_directories("/usr/local/include")
include_directories("/usr/include")
include_directories("/usr/local/include/eigen3")

# include Google Test
set(gtest_SOURCE_DIR "~/src/gtest")
include_directories("${gtest_SOURCE_DIR}/include" ${gtest_SOURCE_DIR})
link_directories("${gtest_SOURCE_DIR}")

set(COMMON_SRCS common/Track.h
    common/Track.cc
    common/TrackMetadata.h)

set(IO_SRCS 
    io/TrackIO.cc
    io/TrackIO.h
    io/Traits.h
    io/BamIO.h
    io/BamIO.cc)
 
set(UTIL_SRCS util/SComm.h
    util/SComm.cc 
    util/SysTime.h
    util/SysTime.cc 
    util/mutil.h
    util/mutil.cc
    util/FastaParser.h)

set(BASE_SRCS base/BEDelement.h
    base/ProgressBar.h 
    base/ProgressBar.cc
    base/WIG.h
    base/WIG.cc
    base/CommandLineParser.h
    base/ErrorHandling.h
    base/ErrorHandling.cc 
    base/EST.h
    base/FileParser.h
    base/FileParser.cc
    base/GFF.h
    base/Interval.h
    base/RandomStuff.h
    base/ReportManager.h
    base/StringUtil.h
    base/StringUtil.cc
    base/SVector.h
    base/DNASequence.h
    base/Common.h
    base/Types.h
    base/NucConversion.h
    base/Hdf5Util.h
    base/Hdf5Util.cc
    base/RefCount.h
    base/Log.h
    base/Log.cc
    base/FileUtil.h
    base/FileUtil.cc)

set(DATA_SRCS data/TrackData.h
    data/ChrData.h
    data/SeqData.h)

set(ANALYSIS_SRCS analysis/Genome.cc 
    analysis/Genome.h
    analysis/Chromosome.h
    analysis/Chromosome.cc
    analysis/GenomeSeq.h
    analysis/GenomeSeq.cc
    analysis/AnalysisBase.h)

SOURCE_GROUP(Util FILES ${UTIL_SRCS})
SOURCE_GROUP(Analysis FILES ${ANALYSIS_SRCS})
SOURCE_GROUP(Data FILES ${DATA_SRCS}) 
SOURCE_GROUP(IO FILES ${IO_SRCS})
SOURCE_GROUP(Base FILES ${BASE_SRCS})
SOURCE_GROUP(Common FILES ${COMMON_SRCS})

# libraries
add_library(util ${COMMON_SRCS} ${IO_SRCS} ${HDF5_SRCS} ${UTIL_SRCS} ${BASE_SRCS} ${DATA_SRCS} ${ANALYSIS_SRCS})
target_link_libraries(util hdf5 bamtools)
# binaries
#add_executable(LoadSeq tool/LoadSeq.cc)

add_executable(LoadData tool/LoadData.cc)
target_link_libraries(LoadData util)

add_executable(LoadSeq tool/LoadSeq.cc)
target_link_libraries(LoadSeq util)

add_executable(RnaSeqCoverage tool/RnaSeqCoverage.cc)
target_link_libraries(RnaSeqCoverage util)

# testing
function(cxx_test name)
  add_executable(${name} "test/${name}.cc")
  target_link_libraries(${name} util gtest)
endfunction()

#cxx_test(GenomeTest)
#cxx_test(ChromosomeTest)
cxx_test(HDFTrackTest)
cxx_test(BamIOTest)
cxx_test(LoadDataTest)
cxx_test(HMMTest)
