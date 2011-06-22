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

#ifndef HDF_REF_GROUP_H_
#define HDF_REF_GROUP_H_

#include "HDFAtom.h"
#include "HDFArray.h"
#include "HDFGroup.h"
#include "datastructures/saf/RefGroup.h"

class HDFRefGroupGroup {
 public:
	HDFGroup refGroup;
	HDFAtom<uint32_t> nRow;
	HDFAtom<string>   masterDatasetAtom;
	HDFArray<uint32_t>  idArray;
	HDFAtom<uint64_t> idArrayLastRow;
	HDFArray<string> pathArray;

	HDFArray<uint32_t> refInfoIdArray;

	~HDFRefGroupGroup() {
		refGroup.Close();
	}
	
	int Initialize(HDFGroup &rootGroup) {
		refGroup.Initialize(rootGroup.group, "RefGroup");
		
		if (masterDatasetAtom.Initialize(refGroup.group, "MasterDataset") == 0) { return 0; }
		if (nRow.Initialize(refGroup, "nRow") == 0) { return 0; }
		if (idArray.Initialize(refGroup, "ID") == 0) { return 0; }
		if (pathArray.Initialize(refGroup, "Path") == 0) { return 0; }
		if (refInfoIdArray.Initialize(refGroup, "RefInfoID") == 0) { return 0; }
		
		return 1;
	}

	void Read(RefGroup &refGroup) {
		masterDatasetAtom.Read(refGroup.masterDataset);
		nRow.Read(refGroup.nRow);
		
		refGroup.path.resize(refGroup.nRow);
		int i;
		int pathNElem = pathArray.arrayLength;
		for (i = 0; i < pathNElem; i++) {
			pathArray.Read(i, i+1, &refGroup.path[i]);
		}

		int idArrayNElem = idArray.arrayLength;
		refGroup.id.resize(idArrayNElem);
		idArray.Read(0, idArrayNElem, &refGroup.id[0]);

		int refIDNElem = refInfoIdArray.arrayLength;
		refGroup.refInfoId.resize(refIDNElem);
		refInfoIdArray.Read(0, refIDNElem, &refGroup.refInfoId[0]);
	}
};

#endif
