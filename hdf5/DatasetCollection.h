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

#ifndef DATA_HDF_DATASET_COLLECTION_H_
#define DATA_HDF_DATASET_COLLECTION_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "hdf5/HDFGroup.h"
#include "hdf5/HDFData.h"


using namespace std;

class DatasetCollection {
 public:
	vector<string> fieldNames;
	map<string,bool> includedFields;
	map<string,bool> requiredFields;

	void MakeFieldRequired(string &fieldName) {
		includedFields[fieldName] = true;
		requiredFields[fieldName] = true;
	}

	void MakeFieldOptional(string &fieldName) {
		includedFields[fieldName] = true;
		requiredFields[fieldName] = false;
	}

	void InitializeAllFields(bool value) {
		int f;
		for (f = 0; f < fieldNames.size(); f++ ) {
			includedFields[fieldNames[f]] = value;
		}
	}

	void InitializeFields(vector<string> &fieldList) {
		int i;
		for (i = 0; i < fieldList.size(); i++) {
			includedFields[fieldList[i]] = true;
		}
	}

	void InitializeFields(vector<char*> &fieldList) {
		int i;
		InitializeAllFields(false);
		for (i = 0; i < fieldList.size(); i++) {
			includedFields[fieldList[i]] = true;
		}
	}

	int IncludeField(string fieldName) {
		if (includedFields.find(fieldName) == includedFields.end()) {
			return 0;
		}	
		else {
			includedFields[fieldName] = true;
		}
	}

	bool ContainsField(string fieldName) {
		int f;
		for (f = 0; f < fieldNames.size(); f++) {
			if (fieldNames[f] == fieldName) return true;
		}
		return false;
	}

	template <typename T_Dataset>
	bool InitializeDataset(HDFGroup &group, T_Dataset &dataset, string datasetName) {
		//
		// Perform initialization of the dataset in a way that keep track
		// of which datasets in the collection are present.
		//
		if (includedFields[datasetName]) {
			if (dataset.Initialize(group, datasetName) == false) {
				if (requiredFields[datasetName]) {
					return false;
				}
				else {
					//
					// This field was supposed to be included but it either does
					// not exist or there was a problem otherwise in creating
					// it.  Don't try and read from it later on.
					//
					includedFields[datasetName] = false;
				}
			}
		}
		return true;
	}

};
	
		

#endif
