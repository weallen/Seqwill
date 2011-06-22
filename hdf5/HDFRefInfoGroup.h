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

#ifndef DATA_HDF_HDF_REF_INFO
#define DATA_HDF_HDF_REF_INFO

#include "datastructures/saf/RefInfo.h"

class HDFRefInfoGroup {
 public:
	HDFGroup refInfoGroup;
	HDFArray<string> fullNameArray;
	HDFArray<uint32_t> idArray;
	HDFArray<uint32_t> lengthArray;
	HDFArray<string> md5Array;
	
	int Initialize(HDFGroup &parentGroup) {
		if (refInfoGroup.Initialize(parentGroup.group, "RefInfo") == 0) {
			return 0;
		}
		
		if (fullNameArray.Initialize(refInfoGroup, "FullName") == 0) { return 0;}
		if (idArray.Initialize(refInfoGroup,"ID") == 0) { return 0;}
		if (lengthArray.Initialize(refInfoGroup, "Length") == 0) { return 0;}
		if (md5Array.Initialize(refInfoGroup, "MD5") == 0) { return 0;}

		return 1;
	}

	~HDFRefInfoGroup() {
		refInfoGroup.Close();
	}

	void Read(RefInfo &refInfo) {

		refInfo.fullName.resize(fullNameArray.arrayLength);
		refInfo.id.resize(idArray.arrayLength);
		refInfo.length.resize(lengthArray.arrayLength);
		refInfo.md5.resize(md5Array.arrayLength);
		if (refInfo.fullName.size() != refInfo.id.size() or
				refInfo.id.size() != refInfo.length.size() or
				refInfo.length.size() != refInfo.md5.size()) {
			cout << "Error with the RefInfo group in a cmp.h5 file.  The datasets " << endl
					 << "are of different lengths but should be the same." << endl;
			exit(1);
		}
		int i;
		for (i = 0; i < refInfo.fullName.size();i++) {
			fullNameArray.Read(i,i+1,&refInfo.fullName[i]);
		}
		for (i = 0; i < refInfo.md5.size(); i++) {
			md5Array.Read(i, i+1, &refInfo.md5[i]);
		}
		lengthArray.Read(0, refInfo.length.size(), &refInfo.length[0]);
		idArray.Read(0, refInfo.id.size(), &refInfo.id[0]);
	}
};


#endif
