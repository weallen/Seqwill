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

#ifndef DATA_HDF_HDF_PULSE_DATA_FILE_H_
#define DATA_HDF_HDF_PULSE_DATA_FILE_H_

#include "data/hdf/HDFGroup.h"
#include "HDFScanDataReader.h"

class HDFPulseDataFile {
 public:
	H5File hdfBasFile;
	HDFGroup pulseDataGroup;
	HDFGroup rootGroup;
	HDFGroup *rootGroupPtr;
	string pulseDataGroupName;
	HDFScanDataReader scanDataReader;
	bool useScanData;
	bool closeFileOnExit;
	int  maxAllocNElements;

	void CheckMemoryAllocation(long allocSize, long allocLimit, const char *fieldName = NULL) {
		if (allocSize > allocLimit) {
			if (fieldName == NULL) {
				cout << "Allocating too large of memory" << endl;
			}
			else {
				cout << "ERROR! Reading the dataset " << fieldName << " will use too much memory." << endl;
				cout << "The pls/bas file is too large, exiting." << endl;
			}
			exit(1);
		}
	}
	
	HDFPulseDataFile() {
		pulseDataGroupName = "PulseData";
		useScanData        = false;
		closeFileOnExit    = false;
		maxAllocNElements  = INT_MAX;
	}

	int OpenHDFFile(string fileName) {
		try {
			H5::FileAccPropList propList;
			hdfBasFile.openFile(fileName.c_str(), H5F_ACC_RDONLY, propList);	
		}
		catch (Exception &e) {
			cout << e.getDetailMsg() << endl;
			cout << "Could not open hdf file" << fileName << ", exiting." << endl;
			exit(0);
			return 0;
		}
		closeFileOnExit = true;
		return 1;
	}

	//
	// All pulse data files contain the "PulseData" group name.
	// 
	//
	int InitializePulseDataFile(string fileName) {
		
		if (OpenHDFFile(fileName) == 0) return 0;
		return 1;
	}

	int Initialize(string fileName) {
		if (InitializePulseDataFile(fileName) == 0) {
			return 0;
		}
		//
		// The pulse group is contained directly below the root group.
		//
		if (rootGroup.Initialize(hdfBasFile, "/") == 0) {
			return 0;
		}
		rootGroupPtr = &rootGroup;
		return Initialize();
	}

	//
	// Initialize inside another open group.
	//
	int Initialize(HDFGroup *rootGroupP) {
		rootGroupPtr = rootGroupP;
		return Initialize();
	}

	//
	// Initialize all fields 
	int Initialize() {
		if (InitializePulseGroup() == 0) return 0;
		if (rootGroupPtr->ContainsObject("ScanData")) {
			if (scanDataReader.Initialize(rootGroupPtr) == 0) {
				return 0;
			}
			else {
				useScanData = true;
			}
		}


		return 1;
	}

	int InitializePulseGroup() {
		if (pulseDataGroup.Initialize(rootGroupPtr->group, pulseDataGroupName) == 0) return 0;
		return 1;
	}
		
	void Close() {
		if (useScanData) {
			scanDataReader.Close();
		}
		
		pulseDataGroup.Close();
		if (rootGroupPtr == &rootGroup) {
			rootGroup.Close();
		}
		/*
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_FILE) << " open files upon closing." <<endl;
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_DATASET) << " open datasets upon closing." <<endl;
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_GROUP) << " open groups upon closing." <<endl;
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_DATATYPE) << " open datatypes upon closing." <<endl;
		cout << "there are " <<  hdfBasFile.getObjCount(H5F_OBJ_ATTR) << " open attributes upon closing." <<endl;
		*/
		if (closeFileOnExit) {
			hdfBasFile.close();
		}


	}

};

#endif
