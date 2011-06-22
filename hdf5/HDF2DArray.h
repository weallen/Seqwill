// Copyright (c) 2010, Pacific Biosciences of California, Inc.
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//     * Neither the name of Pacific Biosciences nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef DATA_HDF_HDF_2DARRAY_H_
#define DATA_HDF_HDF_2DARRAY_H_

#include <H5Cpp.h>
#include "hdf5/HDFConfig.h"
#include "hdf5/HDFData.h"
#include "hdf5/BufferedHDF2DArray.h"
#include <assert.h>
#include <iostream>

using namespace std;
using namespace H5;


/*
 *
 * Implementation of a 2-D array for IO from an HDF array.
 * This is templated, but specialized for a few data types, so that 
 * the HDF data types do not need to be specified by somebody when reading.
 *
 * Currently no support exists for reading non-contiguous blocks of data, and
 * the main intended use is to read in increments of rows.

 int main(int argc, char* argv[]) {
	if (argc < 1) {
		cout << "usage: testHDFReading hdfFile" << endl;
		exit(0);
	}

	string hdfFileName = argv[1];
	
	H5File hdfFile;
	hdfFile.openFile(hdfFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	HDF2DArray<uint16_t> xyArray;
	xyArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/HoleXY");
	int curX = 0;
	xyArray.Read(curX, curX + 1, 0, 2, holeXY);

	or, to read a row:
	xyArray.Read(curX, curX+1, holeXY);

 *
 */
template<typename T>
class HDF2DArray : public BufferedHDF2DArray<T> {
 public:
	void WriteRow(const T*data, int dataLength) {
		this->writeBuffer = (T*) data;
		this->bufferIndex = dataLength;
		this->bufferSize  = dataLength;
		this->Flush();
		//
		// Reset status of buffer so that no methods are tricked into
		// thinking this is a valid pointer.
		//
		this->writeBuffer = NULL;
		this->bufferIndex = 0;
		this->bufferSize  = 0;
	}

};

#endif
