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

#ifndef DATA_HDF_HDF_ATOM_H_
#define DATA_HDF_HDF_ATOM_H_

#include "H5Cpp.h"
#include "data/hdf/HDFConfig.h"
#include "data/hdf/HDFGroup.h"
#include "data/hdf/HDFData.h"
#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>

using namespace std;
using namespace H5;

template<typename T>
class HDFAtom : public HDFData {
 public:
	Attribute attribute;

	bool initialized;
	HDFAtom() {
		initialized = false;
	}
	~HDFAtom() {
		if (initialized) {
			attribute.close();
		}
	}
  int Initialize(H5Object &object, string attributeName) {
		attribute = object.openAttribute(attributeName.c_str());
		initialized = true;
		return 1;
	}
	
	int Initialize(HDFGroup &group, string attributeName) {
		return Initialize(group.group, attributeName);
	}
	
	int Initialize(HDFData &data, string attributeName) {
		return Initialize(data.dataset, attributeName);
	}

  int Initialize(Group &object, string attributeName) {
		attribute = object.openAttribute(attributeName.c_str());
		initialized  = true;
		return 1;
	}

	int Initialize(H5File &hdfFile, string groupName, string attributeName) {
		HDFGroup group;
		group.Initialize(hdfFile, groupName);
		attribute = group.group.openAttribute(attributeName.c_str());
		initialized = true;
		return 1;
	}

	//
	// This handles creation of all non-string types.  A specialization
	// for strings is provided below.
	//
	void Create(H5Object &object, string atomName) {
		hsize_t defaultDims[] = {1};
		DataSpace defaultDataSpace(1, defaultDims);
		TypedCreate(object, atomName, defaultDataSpace);
	}
	

	void Create(H5Object &object, string name, string value) {
		StrType strType(0, value.size());
		hsize_t defaultDims[] = {};
		//		DataSpace defaultDataSpace(1, defaultDims);
		attribute = object.createAttribute(name.c_str(), strType, DataSpace(0,NULL));
		initialized = true;
		attribute.write(strType, value.c_str());
	}
	
	void TypedCreate(H5Object &object, string &atomName, DataSpace &dataSpace) {
		assert("Calling HDFAtom<T>::ypedCreate on an unsupported type" == 0);
	}
	
	void Write(T value) {
		assert("Calling HDFAtom<T>::Write on an unsupported type" == 0);
	}

	void Read(T& value) {
		assert("Calling read on an unsupported type!" == 0);
	}

};

//
// Special create for strings.  Since this uses a StrType for the
// typename rather than specifying a PredType, it mertis its own
// function.
//

template<>
void HDFAtom<string>::Create(H5Object &object, string atomName) {
	StrType strType(0, H5T_VARIABLE);
	hsize_t defaultDims[] = {1};
	DataSpace defaultDataSpace(1, defaultDims);
	attribute = object.createAttribute(atomName.c_str(), strType, defaultDataSpace);
	initialized= true;
}


#define MAKE_TYPED_CREATE(T, predtype) template<> \
	void HDFAtom<T>::TypedCreate(H5Object &object, string &atomName, DataSpace &defaultDataSpace) {				\
  attribute = object.createAttribute(atomName.c_str(), (predtype), defaultDataSpace );	\
}


MAKE_TYPED_CREATE(int, PredType::NATIVE_INT)
MAKE_TYPED_CREATE(unsigned int, PredType::NATIVE_UINT)
MAKE_TYPED_CREATE(unsigned char, PredType::NATIVE_UINT8)
MAKE_TYPED_CREATE(char, PredType::NATIVE_INT8)
MAKE_TYPED_CREATE(float, PredType::NATIVE_FLOAT)
MAKE_TYPED_CREATE(uint64_t, PredType::STD_I64LE)

template<>
void HDFAtom<string>::Write(string value) {
	StrType strType(0, value.size());
	attribute.write(strType, H5std_string(value.c_str()));
}

template<>
void HDFAtom<uint64_t>::Write(uint64_t value) {
	attribute.write( PredType::STD_I64LE, &value);
}

template<>
void HDFAtom<int>::Write(int value) {
	attribute.write( PredType::NATIVE_INT, &value);
}

template<>
void HDFAtom<unsigned int>::Write(unsigned int value) {
	attribute.write( PredType::NATIVE_INT, &value);
}

template<>
void HDFAtom<unsigned char>::Write(unsigned char value) {
	attribute.write( PredType::NATIVE_UINT8, &value);
}

template<>
void HDFAtom<char>::Write(char value) {
	attribute.write( PredType::NATIVE_INT8, &value);
}

template<>
void HDFAtom<float>::Write(float value) {
	attribute.write( PredType::NATIVE_FLOAT, &value);
}


template<>
void HDFAtom<string>::Read(string &value) {
	/*
	 * Read in a string that has been stored either as an array or a
	 * variable length string.  To decide which, query the
	 * isVariableStr() option.
	 */
	StrType stringType = attribute.getStrType();
	bool stringIsVariableLength = stringType.isVariableStr();
	if (stringIsVariableLength) 
		attribute.read(stringType, value);
	else {
		hsize_t stsize = attribute.getStorageSize();
		value.resize(stsize);
		//		char *valueStr = new char[stsize+1];
		attribute.read(stringType, &value[0]);
		if (stsize > 0 and value[stsize-1] == '\0') {
			value.resize(stsize-1);
		}
		//		valueStr[stsize] = '\0';
		//		value = valueStr;
		// This read an extra '\0', which is handled by the string class
		//		if (stsize > 0) {
			//			value = valueStr;
		//			delete[] valueStr;
		//		}
	}
}

template<>
void HDFAtom<int>::Read(int &value) {
	DataType intType(PredType::NATIVE_INT);
	attribute.read(intType, &value);
}

template<>
void HDFAtom<uint64_t>::Read(uint64_t &value) {
	DataType intType(PredType::STD_I64LE);
	attribute.read(intType, &value);
}

template<>
void HDFAtom<unsigned int>::Read(unsigned int &value) {
	DataType uintType(PredType::NATIVE_UINT);
	attribute.read(uintType, &value);
}

template<>
void HDFAtom<float>::Read(float &value) {
	DataType type(PredType::NATIVE_FLOAT);
	attribute.read(type, &value);
}

template<>
void HDFAtom<vector<string> >::Read(vector<string> &values) {
	string value;

	/*
	 * This attribute is an array of strings. They are read in by
	 * storing pointers to strings in memory controlled by HDF.  To read
	 * the strings, read the pointers into a temporary array, then copy
	 * those strings to the values array. This way when the values array
	 * is destroyed, it will not try and get rid of space that is under
	 * HDF control.
	 */
	DataSpace attributeSpace = attribute.getSpace();
	hsize_t nPoints;
	nPoints = attributeSpace.getSelectNpoints();
	DataType attrType = attribute.getDataType(); // necessary for attr.read()
	hsize_t stsize = attribute.getStorageSize();

	// Declare and initialize vector of pointers to string attribute list.
	//	if (nPoints > 1) {
	vector<char*> ptrsToHDFControlledMemory;
	ptrsToHDFControlledMemory.resize(nPoints);
	// Copy the pointers.
	attribute.read(attrType, &ptrsToHDFControlledMemory[0]);
	// Copy the strings into memory the main program has control over.
	unsigned int i;
	for (i = 0; i < ptrsToHDFControlledMemory.size(); i++ ){
		values.push_back(ptrsToHDFControlledMemory[i]);
		free(ptrsToHDFControlledMemory[i]);
	}
		/*
	else if (nPoints == 1) {
		vector<char* > ptrToControlledMemory;
		ptrToControlledMemory.resize(nPoints);
		ptrToControlledMemory[0] = new char[stsize];
		attribute.read(attrType, ptrToControlledMemory[0]);
		values.push_back(ptrToControlledMemory[0]);
		delete ptrToControlledMemory[0];
	}
		*/

}

#endif
