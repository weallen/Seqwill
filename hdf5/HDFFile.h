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

#ifndef DATA_HDF_HDF_FILE_H_
#define DATA_HDF_HDF_FILE_H_

#include <H5Cpp.h>
#include "hdf5/HDFConfig.h"
#include <iostream>
#include <string>
#include <vector>

using namespace H5;
using namespace std;


class HDFFile {
 public:
	// Make this public for easy access.
	H5File *hdfFile;	
	Group rootGroup;
	HDFFile() {
		hdfFile = NULL;
	}

	void Create(string fileName) {
		//		hdfFile.openFile(fileName.c_str(), 
		try {
			FileCreatPropList filePropList;
			hsize_t ub = filePropList.getUserblock();
			filePropList.setUserblock(512);
			hdfFile = new H5File(fileName.c_str(), H5F_ACC_TRUNC, filePropList);
			rootGroup = hdfFile->openGroup("/");
		}
		catch (FileIException fileException) {
			cout << "Error creating file " << fileName << endl;
		}
		return;
	}

	int OpenReadonly(string &fileName) {
		hdfFile->openFile(fileName.c_str(), H5F_ACC_RDONLY);
		return 1;
	}

	/*
		Add a group with an artibrary path name to the file.
		Perhaps this is an abuse of the HDF exception system, but there
		wasn't an easy way to check for existence of a group in the file,
		so simply catching the exceptions that are thrown is one approach.
	*/


	void OpenGroup(string groupName, Group &group) {
		group = hdfFile->openGroup(groupName.c_str());
	}

	void Close() {
		if (hdfFile != NULL) {
			hdfFile->close();
			delete hdfFile;
			hdfFile = NULL;
		}
	}

};


#endif
