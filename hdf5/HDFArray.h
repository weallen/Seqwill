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

#ifndef DATA_HDF_HDF_ARRAY_H_
#define DATA_HDF_HDF_ARRAY_H_

#include <hdf5.h>
#include <H5Cpp.h>
#include "hdf5/HDFConfig.h"
#include "hdf5/HDFData.h"
#include "hdf5/BufferedHDFArray.h"
#include <assert.h>
#include <iostream>
#include <DNASequence.h>
#include <FASTQSequence.h>

using namespace std;
using namespace H5;

/*
 *
 * Implementation of a 1-D array for IO from an HDF array.
 * This is templated, but specialized for a few data types, so that 
 * the HDF data types do not need to be specified by anybody.
 *
 *  Two examples of the usage of this class follow:
 *
 *	HDFArray<int> nElemArray;
 * 	nElemArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/NumEvent");
 *  nElemArray.Read(i, i+1, &nElem);
 *
 * 	HDFArray<unsigned char> qualArray;
 *	qualArray.Initialize(hdfFile, "PulseData/BaseCalls/QualityValue");
 *  qualArray.Read(cur, cur + nElem, qual);
 *
 */

template<typename T>
class HDFArray : public BufferedHDFArray<T> {
 public:

 HDFArray() : BufferedHDFArray<T>() {}
 HDFArray(CommonFG* _container, string _datasetName) : BufferedHDFArray<T>(_container, _datasetName) {}

	/*
	 *  An unbuffered write is simply a write immediately followed by a flush. 
	 */
	void WriteToPos(const T*data, int dataLength, UInt writePos) {
		this->writeBuffer = (T*) data;
		this->bufferIndex = dataLength;
		this->bufferSize  = dataLength;
		this->Flush(false, writePos);
		ResetBuffer();
	}

	void ResetBuffer() {
		this->writeBuffer = NULL;
		this->bufferIndex = 0;
		this->bufferSize  = 0;
	}

	void Write(const T *data, int dataLength) {
		this->writeBuffer = (T*) data;
		this->bufferIndex = dataLength;
		this->bufferSize  = dataLength;
		this->Flush();
		//
		// Reset status of buffer so that no methods are tricked into
		// thinking this is a valid pointer.
		//
		ResetBuffer();
	}
	
	~HDFArray() {} 
	

};


#endif
