cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CC "gcc-4.2")
set(CXX "g++-4.2")

set(CMAKE_C_COMPILER ${CC})

set(CMAKE_CXX_COMPILER ${CXX})
project(Seqwill)

#set(CMAKE_CXX_FLAGS "-O3 -fopenmp -msse -msse2 -msse3 -msse4.2 -m64 -mfpmath=sse")
set(CMAKE_CXX_FLAGS "-g")
#set(cxx_exception_flags "-fexceptions")
#set(cxx_no_exception_flags "-fnoexceptions")

set(EXECUTABLE_OUTPUT_PATH bin)
include("DevSettings.txt")

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
	      base/BEDelement.cc
   base/MatrixUtil.h
   base/MatrixUtil.cc
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

set(ANALYSIS_SRCS 
	  analysis/Nucleosome.h
		analysis/Nucleosome.cc
		analysis/HSMM.cc
		analysis/HSMM.h
		analysis/Genome.cc
    analysis/MedipNormalize.cc
    analysis/MedipNormalize.h
    analysis/gsl_addon.cc
    analysis/gsl_addon.h 
    analysis/Kmeans.h
    analysis/Kmeans.cc
    analysis/HMM.h
    analysis/HMM.cc
    analysis/Dist.h
    analysis/Dist.cc 
    analysis/Random.h
    analysis/Random.cc
    analysis/Genome.h
    analysis/Chromosome.h
    analysis/Chromosome.cc
    analysis/KDE.h
    analysis/KDE.cc
    analysis/AnalysisBase.h)

SOURCE_GROUP(Util FILES ${UTIL_SRCS})
SOURCE_GROUP(Analysis FILES ${ANALYSIS_SRCS})
SOURCE_GROUP(Data FILES ${DATA_SRCS}) 
SOURCE_GROUP(IO FILES ${IO_SRCS})
SOURCE_GROUP(Base FILES ${BASE_SRCS})
SOURCE_GROUP(Common FILES ${COMMON_SRCS})

set(ALL_SRCS ${COMMON_SRCS} ${IO_SRCS} ${HDF5_SRCS} ${UTIL_SRCS} ${BASE_SRCS})

# libraries
add_library(SWLib ${COMMON_SRCS} ${IO_SRCS} ${HDF5_SRCS} ${UTIL_SRCS} ${BASE_SRCS} ${DATA_SRCS} ${ANALYSIS_SRCS})
target_link_libraries(SWLib hdf5 bamtools gsl gslcblas tbb hoard)
# binaries

function(tool name) 
  add_executable(${name} tool/${name}.cc)
  target_link_libraries(${name} SWLib)
endfunction()

function(mpi_tool name)
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpicxx)		
  add_executable(${name} tool/${name}.cc)
  target_link_libraries(${name} SWLib boost_mpi boost_serialization mpi)
  set(CMAKE_C_COMPILER ${CC})
  set(CMAKE_CXX_COMPILE ${CXX})
endfunction()

tool(LoadData)

tool(LoadSeq)

tool(GenerateHMMTestData)

tool(ReadStartPileup)

tool(Phasogram)

tool(TrackToWig)

tool(SubTrackToWig)

tool(CountCpGs)

tool(DipNormalize)

tool(DipNormalizeTBB)
#tool(DipNormalizeOMP)
#mpi_tool(DipNormalizeParallel)

tool(CombineTrackFiles)

tool(RPMNormalize)

tool(HMMCompare)

tool(CpGAverage)

tool(NucMap)

tool(NucPosition)

# testing
function(cxx_test name)
  add_executable(${name} "test/${name}.cc")
  target_link_libraries(${name} gtest SWLib)
endfunction()

#cxx_test(GenomeTest)
#cxx_test(ChromosomeTest)
cxx_test(HDFTrackTest)
cxx_test(LoadDataTest)
cxx_test(HMMTest)
cxx_test(DistTest)
cxx_test(BamIOTest)
cxx_test(MedipNormTest)
cxx_test(NucMapTest)
#cxx_test(KmeansTest)
