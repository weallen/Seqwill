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

#ifndef DATA_HDF_HDF_DATA_H_
#define DATA_HDF_HDF_DATA_H_


#include <H5Cpp.h>
#include "hdf5/HDFConfig.h"
#include "hdf5/HDFAttributable.h"
#include <assert.h>
#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace H5;

class HDFData : public HDFAttributable {
 public:
	DataSet   dataset;
	DataSpace dataspace;
	DataSpace sourceSpace, destSpace;
	DataSpace fullSourceSpace;
	bool      fileDataSpaceInitialized;
	CommonFG  *container;
	string    datasetName;
	bool      isInitialized;
	
 HDFData(CommonFG* _container, string _datasetName) {
		HDFData();
		container   = _container;
		datasetName = _datasetName;
	}

	HDFData() {
		fileDataSpaceInitialized = false;
		isInitialized = false;
	}

	bool IsInitialized() {
		return isInitialized;
	}

	int InitializeDataset(CommonFG &hdfFile, string _datasetName) {
		try {
			datasetName = _datasetName;
			dataset   = hdfFile.openDataSet(_datasetName.c_str());
			//			size_t datasetSize = dataset.getStorageSize();
			isInitialized = true;
			fileDataSpaceInitialized = true;
		}
		catch(FileIException &e) {
			cerr << e.getDetailMsg() <<endl;
			return 0;
		}
		catch(GroupIException &e) {
			cerr << e.getDetailMsg() << endl;
			return 0;
		}
		catch(H5::Exception &e) {
			cerr << e.getDetailMsg() << endl;
			return 0;
		}
		StoreAttributeNames(dataset);
		return 1;
	}
	
	void Close() {
		dataspace.close();
		dataset.close();
	}
};


			

#endif
