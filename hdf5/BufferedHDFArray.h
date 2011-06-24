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

#ifndef DATA_HDF_BUFFERED_HDF_ARRAY_H_
#define DATA_HDF_BUFFERED_HDF_ARRAY_H_

#include <hdf5.h>
#include <H5Cpp.h>
#include "hdf5/HDFConfig.h"
#include "hdf5/HDFData.h"
#include "hdf5/HDFGroup.h"
#include "hdf5/HDFWriteBuffer.h"
#include "base/DNASequence.h"
#include <assert.h>
#include <iostream>

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
class BufferedHDFArray : public HDFData, public HDFWriteBuffer<T> {
 protected:
 public:
	hsize_t   nDims;
	hsize_t   *dimSize;
	int       maxDims;
	int       arrayLength;
	/*
	 * Constructor meant to be used for data that will be written.  
	 * This allocates the write buffer.
	 */
 BufferedHDFArray(int pBufferSize=0) : HDFData() {
		nDims = 0;
		dimSize = NULL;
		this->bufferIndex = 0;
		this->InitializeBuffer(pBufferSize);
	}
	

	BufferedHDFArray(CommonFG* _container, string _datasetName) : HDFData(_container, _datasetName) {
		//		HDFArray();
	}

	~BufferedHDFArray() {

		//
		// Clean up the write buffer.
		//
		//		Flush();
		if (dimSize != NULL) {
			delete[] dimSize;
		}
		this->HDFWriteBuffer<T>::~HDFWriteBuffer();
	}

	void Initialize(CommonFG* _container, string _datasetName, int pBufferSize = 0) {
		container   = _container;
		datasetName = _datasetName;
		this->InitializeBuffer(pBufferSize);
	}


	void Write(const T *data, int dataLength) {
		// Fill the buffer with data. When there is overflow, write
		// that out to disk.
		//
		int dataIndex = 0;
		int bufferCapacity;
		int bufferFillSize = 0;
		bool flushBuffer;
		while(dataIndex < dataLength) {
			bufferCapacity = this->bufferSize - this->bufferIndex;
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
	/*	void Write(string &str) {
    StrType tid(0, H5T_FIXED);
		dataset.write(str, tid);
		}*/
	void Flush(bool append=true, UInt writePos = 0) {
		//
		// Flush contents of current buffer to the file.
		//
		if (this->WriteBufferEmpty()) {
			// 
			// There is no data in the buffer, so nothing can be written.
			// HDF does not support empty arrays (as far as I can tell), so
			// don't even bother trying to create the dataspace.
			// 
			return;
		}

		// fetch the current size of the dataspace
		if (fileDataSpaceInitialized == false) {
			Create(this->writeBuffer, writePos + this->bufferIndex, writePos);
			fileDataSpaceInitialized = true;
		}
		else {
			DataSpace fileSpace;
			fileSpace = dataset.getSpace();

			//
			// Load the current size of the array on disk.
			//
			hsize_t fileArraySize[1], blockStart;
			fileArraySize[0] = fileSpace.getSimpleExtentNpoints();
			if (append) {
				blockStart = fileSpace.getSimpleExtentNpoints();
				fileArraySize[0] += this->bufferIndex;
				//
				// Make room in the file for the array.
				//
				dataset.extend(fileArraySize);
			}
			else {
				blockStart = writePos;
				if (blockStart + this->bufferIndex > fileArraySize[0]) {
					fileArraySize[0] = blockStart + this->bufferIndex;
					dataset.extend(fileArraySize);
				}
			}
		
			DataSpace extendedSpace = dataset.getSpace();
			int extendedSize        = extendedSpace.getSimpleExtentNpoints();
			//
			// Configure the proper addressing to append to the array.
			//
			hsize_t dataSize[1];
			hsize_t offset[1];
			dataSize[0] = this->bufferIndex;
			offset[0]   = blockStart;
			extendedSpace.selectHyperslab(H5S_SELECT_SET, dataSize, offset);
			DataSpace memorySpace(1, dataSize);

			//
			// Finally, write out the data.  
			// This uses a generic function which is specialized with
			// templates later on to t
			// memorySpace addresses the entire array in linear format
			// fileSpace addresses the last dataLength blocks of dataset.
			//
			TypedWrite(this->writeBuffer, memorySpace, extendedSpace);
			memorySpace.close();
			extendedSpace.close();
			fileSpace.close();
		}
		// Clear the buffer.
		this->ResetWriteBuffer();
	}

	void TypedWrite(const char **data, const DataSpace &memorySpace, const DataSpace &extendedSpace) {
		StrType varStrType(0,H5T_VARIABLE);
		dataset.write(data, varStrType, memorySpace, extendedSpace);
	}


	void TypedWrite(const T*data, const DataSpace &memorySpace, const DataSpace &extendedSpace) {
		assert("Calling TypedWrite on an unsupported type" == 0);
	}

	void Create(const T *data, int dataLength, int writePos=-1) {
		hsize_t dataSize[]    = {dataLength};
		hsize_t maxDataSize[] = {H5S_UNLIMITED};
		DataSpace fileSpace(1, dataSize, maxDataSize);
		DSetCreatPropList cparms;

		/*
		 * For some reason, chunking must be enabled when creating a dataset
		 * that  has an unlimited dimension.  Of course, this is not
		 * mentioned in the hdf5 c++ documentation, because that
		 * docuemntation was written for people who enjoy learning how to
		 * use an API by reading comments in source code.
		 */
		hsize_t      chunk_dims[1] ={16384};
		cparms.setChunk( 1, chunk_dims );
		TypedCreate(fileSpace, cparms);
		
		//
		// Since TypedCreate created an assigned a dataset, this array is
		// now initialized.  Do the bookkeeping here.
		//
		isInitialized = true;

		if (writePos == -1) {
			// append data
			TypedWrite(data, DataSpace::ALL, DataSpace::ALL);
		}
		else {
			DataSpace extendedSpace = dataset.getSpace();
			int extendedSize        = extendedSpace.getSimpleExtentNpoints();
			//
			// Configure the proper addressing to append to the array.
			//
			hsize_t dataSize[1];
			hsize_t offset[1];
			dataSize[0] = this->bufferIndex;
			offset[0]   = writePos;
			extendedSpace.selectHyperslab(H5S_SELECT_SET, dataSize, offset);
			DataSpace memorySpace(1, dataSize);
			TypedWrite(data, memorySpace, extendedSpace);
			extendedSpace.close();
		}
		fileSpace.close();
	}

	void TypedCreate(DataSpace &fileSpace, DSetCreatPropList &cparms) {
		cout << "DEFAULT typed create " << endl;
		//		dataset = container->createDataSet(datasetName.c_str(), PredType::NATIVE_INT, fileSpace, cparms);
	}

	/*
	 * Initialize for reading.
	 *
	 * Open a dataset in an hdf file. Only call this on datasets that
	 * exist, since this currently handles errors with opening datasets
	 * by ungracefully exiting the program. 
	 */

	int Initialize(HDFGroup &parentGroup, string datasetName) {

		//
		// Make sure writing will not work.
		// 
		this->writeBuffer = NULL;
		this->bufferIndex = 0;
		
		//
		// It's possible that the group may be asked to initialize this
		// dataset when the dataset does not exist.  Check that here.
		//
		if (parentGroup.ContainsObject(datasetName) == 0) {
			return 0;
		}
		
		if (HDFData::InitializeDataset(parentGroup.group, datasetName) == 0) {
			return 0;
		}
		try {
			dataspace = dataset.getSpace();
		}
		catch(H5::DataSetIException &e) { 
			return 0;
		}
		maxDims = MAX_DIMS;
		try {
			nDims = dataspace.getSimpleExtentNdims();
			/*
			 * Prevent abuse of this class for multidimensional IO.
			 */
			if (nDims != 1) {
				cout << "ERROR in HDF format: dataset: " << datasetName << " should be 1-D, but it is not." << endl;
				exit(1);
			}

			/*
			 * Load in the size of this dataset, and make a map to the whole thing.
			 */
			dimSize = new hsize_t[nDims];
			dataspace.getSimpleExtentDims(dimSize);
			if (dimSize[0] == 0) {
				return 0;
			}
			arrayLength = dimSize[0];
			fullSourceSpace = DataSpace(1, dimSize);
			dataspace.close();
		}
		catch(Exception &e) {
			return 0;
		}
		return 1;
	}

	int size() {
		dataspace = dataset.getSpace();
		hsize_t dimSizeArray[1];
		dataspace.getSimpleExtentDims(dimSizeArray);
		dataspace.close();
		return dimSizeArray[0];
	}

	/*
	 * Unspecialized form of read.
	 * Read cannot be called on a type T* that does not have a
	 * specialized template definition.  This is all determined at
	 * compile time.  To ensure this, the following
	 * default definition is provided that gives a nasty warning and
	 * exits the code.
	 */


	void Read(int start, int end, T* dest) {
		assert("ERROR, calling Read with an unsupported type. Use Read(start,end,datatype, dest) instead." == 0);
		exit(0); // this is in case the assert statement is removed.
	}

	/*
	 * Read in type T from the opened dataset from the interval (start,
	 * end].
	 */
	
	void ReadDataset(vector<T> &dest) {
		assert("ERROR, calling ReadDataset with an unsupported type.");
		exit(0); // this is in case the assert statement is removed.
	}


	void Read(int start, int end, DataType typeID, T *dest) {
		if (end - start == 0) {
			return;
		}
		hsize_t memSpaceSize[] = {0};
		memSpaceSize[0] = end - start;
		hsize_t sourceSpaceOffset[] = {0};
		sourceSpaceOffset[0] = start;
		DataSpace destSpace(1, memSpaceSize);
		fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
		dataset.read(dest, typeID, destSpace, fullSourceSpace);
		destSpace.close();
	}

	void ReadCharArray(int start, int end, string*dest) {
		hsize_t memSpaceSize[] = {0};
		memSpaceSize[0] = end - start;
		hsize_t sourceSpaceOffset[] = {0};
		sourceSpaceOffset[0] = start;
		DataSpace destSpace(1, memSpaceSize);
		StrType strType(0, H5T_VARIABLE);
		fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
		vector<char*> tmpStringArray;
		tmpStringArray.resize(end-start);
		dataset.read(&tmpStringArray[0], strType, destSpace, fullSourceSpace);
		int i;
		for (i = 0; i < tmpStringArray.size(); i++) {
			dest[i] = tmpStringArray[i];
		}
		destSpace.close();
	}
		
};

/*
 * Type specializations for some standard types. Use the macro for
 * vanilla specializations (that only require the HDF type ID to be
 * specified). 
 */
#define DEFINE_TYPED_READ_ARRAY(T, Pred) template<>  \
   	void BufferedHDFArray<T>::Read(int start, int end, T* dest) { \
   	Read(start,end, Pred, dest); \
	}


DEFINE_TYPED_READ_ARRAY(int, PredType::NATIVE_INT);
DEFINE_TYPED_READ_ARRAY(char, PredType::NATIVE_INT8);
DEFINE_TYPED_READ_ARRAY(unsigned char, PredType::NATIVE_UINT8);
DEFINE_TYPED_READ_ARRAY(unsigned int, PredType::NATIVE_UINT);
DEFINE_TYPED_READ_ARRAY(uint16_t, PredType::NATIVE_UINT16);
DEFINE_TYPED_READ_ARRAY(float, PredType::NATIVE_FLOAT);
DEFINE_TYPED_READ_ARRAY(char*, PredType::C_S1);


#define DEFINE_TYPED_READ_DATASET(T, Pred) template<>  \
	void BufferedHDFArray<T>::ReadDataset(vector<T>  &dest) {	 \
	dest.resize(arrayLength); \
  Read(0,arrayLength, Pred, &dest[0]);										\
}

DEFINE_TYPED_READ_DATASET(int, PredType::NATIVE_INT);
DEFINE_TYPED_READ_DATASET(char, PredType::NATIVE_INT8);
DEFINE_TYPED_READ_DATASET(unsigned char, PredType::NATIVE_UINT8);
DEFINE_TYPED_READ_DATASET(unsigned int, PredType::NATIVE_UINT);
DEFINE_TYPED_READ_DATASET(uint16_t, PredType::NATIVE_UINT16);
DEFINE_TYPED_READ_DATASET(float, PredType::NATIVE_FLOAT);
DEFINE_TYPED_READ_DATASET(char*, PredType::C_S1);

template<>
void BufferedHDFArray<string>::Read(int start, int end, string*dest) {
	vector<char*> tmpDestCharPtrs;
	if (end == start) return;
	assert(end > start);
	tmpDestCharPtrs.resize(end-start);
	//	ReadCharArray(start, end, (char**) &tmpDestCharPtrs[0]);
	ReadCharArray(start, end, dest);
	//	unsigned int ptrIndex;
	//	for (ptrIndex = 0; ptrIndex < end - start; ptrIndex++) {
	//		dest[ptrIndex] = tmpDestCharPtrs[ptrIndex];
	//	}
	// Unset the values of the tmp so that they are not destructed when
	// removed from the stack.
	//	tmpDestCharPtrs.resize(0);
}

#define DEFINE_TYPED_CREATE_ARRAY(T,Pred) template<> \
	void BufferedHDFArray<T>::TypedCreate(DataSpace &fileSpace,  DSetCreatPropList &cparms) { \
	T zero; zero = 0;\
	cparms.setFillValue(Pred,&zero);\
	dataset = container->createDataSet(datasetName.c_str(), Pred, fileSpace, cparms); \
}
	//	isInitialized = true;																							\

DEFINE_TYPED_CREATE_ARRAY(int, PredType::NATIVE_INT);
DEFINE_TYPED_CREATE_ARRAY(char, PredType::NATIVE_INT8);
DEFINE_TYPED_CREATE_ARRAY(char*, StrType(0,H5T_VARIABLE));
//DEFINE_TYPED_CREATE_ARRAY(string, StrType(0,H5T_VARIABLE));
DEFINE_TYPED_CREATE_ARRAY(unsigned char, PredType::NATIVE_UINT8);
DEFINE_TYPED_CREATE_ARRAY(unsigned int, PredType::NATIVE_UINT);
DEFINE_TYPED_CREATE_ARRAY(float, PredType::NATIVE_FLOAT);
DEFINE_TYPED_CREATE_ARRAY(uint16_t, PredType::NATIVE_UINT16);

template<>
void BufferedHDFArray<string>::TypedCreate(DataSpace &space, DSetCreatPropList &cparms) {
	StrType varStrType(0,H5T_VARIABLE);
	dataset = container->createDataSet(datasetName.c_str(), varStrType, space, cparms);
}

#define DEFINE_TYPED_WRITE_ARRAY(T, Pred) template<>													\
	void BufferedHDFArray<T>::TypedWrite(const T *data, const DataSpace &memorySpace, const DataSpace &fileSpace) {	\
		dataset.write(data, Pred, memorySpace, fileSpace);									\
	}


DEFINE_TYPED_WRITE_ARRAY(int, PredType::NATIVE_INT)
DEFINE_TYPED_WRITE_ARRAY(unsigned int, PredType::NATIVE_UINT)
DEFINE_TYPED_WRITE_ARRAY(unsigned char, PredType::NATIVE_UINT8)
DEFINE_TYPED_WRITE_ARRAY(char, PredType::NATIVE_INT8)
DEFINE_TYPED_WRITE_ARRAY(float, PredType::NATIVE_FLOAT)
DEFINE_TYPED_WRITE_ARRAY(uint16_t, PredType::NATIVE_UINT16)
DEFINE_TYPED_WRITE_ARRAY(char*, StrType(0,H5T_VARIABLE))
DEFINE_TYPED_WRITE_ARRAY(string, StrType(0,H5T_VARIABLE))

/*
 * This is a nonstandard definition because it requires the creation
 * of a special datatype for variable length string type.
 */


#endif
