cmake_minimum_required(VERSION 2.8)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#set(CMAKE_C_COMPILER "clang")

#set(CMAKE_CXX_COMPILER "clang++")
project(Seqwill)

#set(CMAKE_CXX_FLAGS "-Wall -Wshadow")
#set(cxx_exception_flags "-fexceptions")
#set(cxx_no_exception_flags "-fnoexceptions")

set(EXECUTABLE_OUTPUT_PATH bin)

include_directories(".")
include_directories("/usr/local/include/bamtools")
include_directories("/usr/local/include")
include_directories("/usr/include")

# include Google Test
set(gtest_SOURCE_DIR "~/src/googletest-read-only")
include_directories("${gtest_SOURCE_DIR}/include" ${gtest_SOURCE_DIR})
link_directories("${gtest_SOURCE_DIR}")

set(HDF5_SRCS hdf5/BufferedHDF2DArray.h
        hdf5/BufferedHDFArray.h
        hdf5/HDFArray.h
        hdf5/HDFAtom.h
        hdf5/HDFAttributable.h
        hdf5/HDFData.h
        hdf5/HDFFile.h
        hdf5/HDFGroup.h
        hdf5/HDFCommonFG.h
        hdf5/HDFConfig.h
        hdf5/HDFWriteBuffer.h

        hdf5/HDFGenomeFile.h
        hdf5/HDFGenomeTrackData.h
        hdf5/HDFGenomeReader.h
        hdf5/HDFGenomeWriter.h
        hdf5/HDFGenomeChrGroup.h
        hdf5/HDFGenomeRootGroup.h

        hdf5/HDFBasReader.h
        hdf5/HDFBasWriter.h
        hdf5/HDFCCSReader.h

        hdf5/HDFCmpData.h
        hdf5/HDFCmpExperimentGroup.h
        hdf5/HDFCmpFile.h
        hdf5/HDFCmpReader.h
        hdf5/HDFCmpRefAlignmentGroup.h
        hdf5/HDFCmpRootGroup.h
        hdf5/HDFCmpWriter.h

        hdf5/HDFMovieInfoGroup.h

        hdf5/HDFPlsReader.h
        hdf5/HDFPlsWriter.h
        hdf5/HDFPulseDataFile.h

        hdf5/HDFRefGroupGroup.h
        hdf5/HDFRefInfoGroup.h
        hdf5/HDFRegionTableReader.h
        hdf5/HDFRegionTableWriter.h

        hdf5/HDFScanDataReader.h
        hdf5/HDFSMRTSequenceReader.h
        hdf5/HDFZMWReader.h
        hdf5/PlatformId.h)

# libraries
add_library(util ${HDF5_SRCS} util/SComm.cc util/SysTime.cc util/mutil.cc base/ErrorHandling.cc base/FileParser.cc base/StringUtil.cc)
add_library(analysis analysis/Genome.cc analysis/Chromosome.cc analysis/Bam.cc)
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
