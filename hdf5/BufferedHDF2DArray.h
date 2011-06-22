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

#ifndef DATA_BUFFERED_HDF_HDF_2DARRAY_H_
#define DATA_BUFFERED_HDF_HDF_2DARRAY_H_

#include <H5Cpp.h>
#include "hdf5/HDFConfig.h"
#include "hdf5/HDFData.h"
#include "hdf5/HDFGroup.h"
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
	BufferedHDF2DArray<uint16_t> xyArray;
	xyArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/HoleXY");
	int curX = 0;
	xyArray.Read(curX, curX + 1, 0, 2, holeXY);

	or, to read a row:
	xyArray.Read(curX, curX+1, holeXY);

 *
 */
template<typename T>
class BufferedHDF2DArray : public HDFData, public HDFWriteBuffer<T> {
	hsize_t   nDims;
	hsize_t   *dimSize;
	int       maxDims;
	int       rowLength, colLength;
 public:
	unsigned int GetNRows() {
		return colLength;
	}
	unsigned int GetNCols() {
		return rowLength;
	}

 BufferedHDF2DArray(CommonFG *_container, string _datasetName) : HDFData(_container, _datasetName) {
	}

 BufferedHDF2DArray() : HDFData() {
		nDims = 2;
		dimSize =NULL;
		rowLength = -1;
		colLength = -1;
	}

	~BufferedHDF2DArray() {

		//
		// Clean up the write buffer.
		//
		//		Flush();
		if (dimSize != NULL) {
			delete[] dimSize;
		}
		this->HDFWriteBuffer<T>::~HDFWriteBuffer();
	}
	
	int Initialize(CommonFG *_container, string _datasetName, int pRowLength, int pBufferSize=0) {
		//		HDF2DArray(_container, _datasetName);
		container   = _container;
		datasetName = _datasetName;
		nDims       = 2;
		dimSize     = NULL;
		rowLength   = pRowLength;
		this->InitializeBuffer(pBufferSize);
		return 1;
	}

	
	/*
	 * Initialize HDF2D for reading.  No write buffer initialization is
	 * required.  The assumption is that the dataspace is in two
	 * dimensions, and this exits without grace if it is not. 
	 */
	int Initialize(HDFGroup &group, string datasetName) {
		return Initialize(group.group, datasetName);
	}

	int Initialize(CommonFG &commonfg, string datasetName) {
		InitializeDataset(commonfg, datasetName);
		try {
			dataspace = dataset.getSpace();
		}
		catch(H5::DataSetIException &e) { 
			cout << e.getDetailMsg() << endl;
			exit(1);
		}
		maxDims   = MAX_DIMS;
		try {
			nDims     = dataspace.getSimpleExtentNdims();
			/*
			 * Prevent abuse of this class for multidimensional IO.
			 */
			if (nDims != 2) {
				cout << "ERROR in HDF format: dataset: " << datasetName << " should be 1-D, but it is not." << endl;
				exit(1);
			}

			/*
			 * Load in the size of this dataset, and make a map to the whole thing.
			 */
			dimSize = new hsize_t[nDims];
			dataspace.getSimpleExtentDims(dimSize);
			fullSourceSpace = DataSpace(2, dimSize);
			colLength = dimSize[0];
			rowLength = dimSize[1];
		}
		catch(Exception &e) {
			cout << e.getDetailMsg() << endl;
			exit(1);
		}
		return 1;
	}

	int size() {
		assert(nDims == 1);
		dataspace.getSimpleExtentDims(dimSize);
		return dimSize[0];
	}

	/*
	 * Read rows in the range (startX, endX] in to dest.
	 */

	void Read(int startX, int endX, DataType typeID, T*dest) {
		Read(startX, endX, 0, dimSize[1], typeID, dest);
	}
	
	void Read(int startX, int endX, T*dest) {
		Read(startX, endX, 0, dimSize[1], dest);
	}
	/*
	 * This is the non-specialized definition.  Since this should only
	 * operate on specialized types, report an error and bail.
	 */
	void Read(int startX, int endX, int startY, int endY, T* dest) {
		assert("ERROR, calling Read with an unsupported type. Use Read(startx,endx, starty,endy,datatype, dest) instead." == 0);
		exit(0);
	}

	void Read(int startX, int endX, int startY, int endY, DataType typeID, T *dest) {
		hsize_t memSpaceSize[2] = {0, 0};
		memSpaceSize[0] = endX - startX;
		memSpaceSize[1] = endY - startY;
		hsize_t sourceSpaceOffset[2] = {0, 0};
		sourceSpaceOffset[0] = startX;
		sourceSpaceOffset[1] = startY;
		
		DataSpace destSpace(2, memSpaceSize);		
		fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
		dataset.read(dest, typeID, destSpace, fullSourceSpace);
		destSpace.close();
	}
	
	void Create(const T *data, int nRows=1) {
		//
		// Make life easy if the buffer is too small to fit a row --
		// resize it so that rows may be copied and written out in an
		// atomic unit.
		//
		if (this->bufferSize < rowLength) {
			// When the buffer size is greater than 0, the write buffer
			// should exist.
			assert(this->writeBuffer != NULL);
			delete[] this->writeBuffer;
			this->writeBuffer = new T[rowLength];
			this->bufferSize = rowLength;
		}

		hsize_t dataSize[2]    = {nRows, rowLength};
		hsize_t maxDataSize[2] = {H5S_UNLIMITED, rowLength};
		DataSpace fileSpace(2, dataSize, maxDataSize);
		DSetCreatPropList cparms;

		/*
		 * For some reason, chunking must be enabled when creating a dataset
		 * that  has an unlimited dimension.  Of course, this is not
		 * mentioned in the hdf5 c++ documentation, because that
		 * docuemntation was written for people who enjoy learning how to
		 * use an API by reading comments in source code.
		 */
		hsize_t chunkDims[2] ={16384, rowLength};
		cparms.setChunk( 2, chunkDims );
		TypedCreate(fileSpace, cparms);
		TypedWriteRow(this->writeBuffer, DataSpace::ALL, DataSpace::ALL);
		fileSpace.close();
	}

	void TypedCreate(DataSpace &fileSpace, DSetCreatPropList &cparms) {
		assert("Error, calling HDF2DArray<T>::TypedCreate on an unsupported type.  A specialization must be written in HDF2DArray.h" == 0);
	}
	// Append
	void TypedWriteRow(const T*, const DataSpace &memoryDataSpace, const DataSpace &fileDataSpace) {
		assert("Error, calling HDF2DArray<T>::TypedWriteRow on an unsupported type.  A specialization must be written in HDF2DArray.h" == 0);
	}


	/*
	 * This code is copied directly form BufferedHDFArray.  I'm not sure
	 * how to set up the objects nicely to share the code between the
	 * two since the Flush() function is different.  There probably is a
	 * design pattern or simply better way to engineer this, but for now
	 * it's 15 lines of code.
	 */
	 
	void WriteRow(const T *data, int dataLength) {
		// Fill the buffer with data. When there is overflow, write
		// that out to disk.
		//
		int dataIndex = 0;
		int bufferCapacity;
		int bufferFillSize = 0;
		bool flushBuffer;
		while(dataIndex < dataLength) {
			//
			// Compute the capacity of this buffer to fit an integral number
			// of rows into it.
			//
			bufferCapacity = (this->bufferSize / rowLength)*rowLength - this->bufferIndex;
			flushBuffer = false;
			if (bufferCapacity  > dataLength - dataIndex) {
				bufferFillSize = dataLength - dataIndex;
			}
			else {
				bufferFillSize = bufferCapacity;
				flushBuffer = true;
			}
			memcpy((void*) &this->writeBuffer[this->bufferIndex], (void*) &data[dataIndex], sizeof(T)*bufferFillSize);
			dataIndex   += bufferFillSize;
			this->bufferIndex += bufferFillSize;
			if (flushBuffer) {
				Flush();
			}
		}
	}

	void Flush() {
		int numRows = this->bufferIndex/ rowLength;
		if (fileDataSpaceInitialized == false) {
			if (numRows > 0) {
				Create(this->writeBuffer, numRows);
				fileDataSpaceInitialized = true;
			}
 		}
		else {
			if (numRows > 0) {
				DataSpace fileSpace;
				fileSpace = dataset.getSpace();
			
				//
				// Load the current size of the array on disk.
				//
				hsize_t fileArraySize[2], fileArrayMaxSize[2], blockStart[2];
				fileSpace.getSimpleExtentDims(fileArraySize, fileArrayMaxSize);
				// Save this for later to determine the offsets
				blockStart[0] = fileArraySize[0];
				blockStart[1] = fileArraySize[1];

				// increment the size by one row.
				fileArraySize[0]+= numRows;

				//
				// Make room in the file for the array.
				//
				dataset.extend(fileArraySize);
		
				DataSpace extendedSpace = dataset.getSpace();
				//
				// Store the newly dimensioned dataspaces.
				//
				fileSpace.getSimpleExtentDims(fileArraySize, fileArrayMaxSize);			
				int extendedSize = extendedSpace.getSimpleExtentNpoints();
				//
				// Configure the proper addressing to append to the array.
				//
				hsize_t dataSize[2];
				dataSize[0] = numRows;
				dataSize[1] = rowLength;
				hsize_t offset[2];
				offset[0] = blockStart[0];
				offset[1] = 0;
				extendedSpace.selectHyperslab(H5S_SELECT_SET, dataSize, offset);
				DataSpace memorySpace(2, dataSize);

				//
				// Finally, write out the data.  
				// This uses a generic function which is specialized with
				// templates later on to t
				// memorySpace addresses the entire array in linear format
				// fileSpace addresses the last dataLength blocks of dataset.
				//
				TypedWriteRow(this->writeBuffer, memorySpace, extendedSpace);
				memorySpace.close();
				extendedSpace.close();
				fileSpace.close();
			}
		}
		this->ResetWriteBuffer();
	}

};

#define DEFINE_TYPED_WRITE_ROW(T, Pred) template<>\
void BufferedHDF2DArray<T>::TypedWriteRow(const T *data, const DataSpace &memorySpace, const DataSpace &fileSpace) {\
	dataset.write(data, Pred, memorySpace, fileSpace);\
}

DEFINE_TYPED_WRITE_ROW(int, PredType::NATIVE_INT);
DEFINE_TYPED_WRITE_ROW(unsigned int, PredType::NATIVE_UINT);
DEFINE_TYPED_WRITE_ROW(unsigned char, PredType::NATIVE_UINT8);
DEFINE_TYPED_WRITE_ROW(uint16_t, PredType::NATIVE_UINT16);
DEFINE_TYPED_WRITE_ROW(int16_t, PredType::NATIVE_INT16);


#define DEFINE_READ(T, Pred) template<>\
void BufferedHDF2DArray<T>::Read(int startX, int endX, int startY, int endY, T* dest) {\
	Read(startX, endX, startY, endY, Pred, dest);\
}


DEFINE_READ(int, PredType::NATIVE_INT);
DEFINE_READ(unsigned int, PredType::NATIVE_UINT);
DEFINE_READ(char, PredType::NATIVE_INT8);
DEFINE_READ(unsigned char, PredType::NATIVE_UINT8);
DEFINE_READ(uint16_t, PredType::NATIVE_UINT16);
DEFINE_READ(int16_t, PredType::NATIVE_INT16);

#define DEFINE_TYPED_CREATE(T, Pred)template<>\
void BufferedHDF2DArray<T>::TypedCreate(DataSpace &fileSpace, DSetCreatPropList &cparms) {\
	dataset = container->createDataSet(datasetName.c_str(), Pred, fileSpace, cparms);\
}

DEFINE_TYPED_CREATE(int, PredType::NATIVE_INT)
DEFINE_TYPED_CREATE(unsigned int, PredType::NATIVE_UINT)
DEFINE_TYPED_CREATE(char, PredType::NATIVE_INT8)
DEFINE_TYPED_CREATE(unsigned char, PredType::NATIVE_UINT8)
DEFINE_TYPED_CREATE(uint16_t, PredType::NATIVE_UINT16)
DEFINE_TYPED_CREATE(int16_t, PredType::NATIVE_INT16)


#endif
