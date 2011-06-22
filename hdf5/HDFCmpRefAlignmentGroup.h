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

#ifndef DATA_HDF_HDF_CMP_REF_ALIGNMENT_GROUP_H_
#define DATA_HDF_HDF_CMP_REF_ALIGNMENT_GROUP_H_

#include "H5Cpp.h"
#include "data/hdf/HDFData.h"
#include "data/hdf/HDFGroup.h"
#include "data/hdf/HDFCmpExperimentGroup.h"
#include <map>
#include <string>
using namespace std;

class HDFCmpRefAlignmentGroup {
 public:
	HDFGroup  refGroup;
	string refGroupName;
	vector<HDFCmpExperimentGroup*> readGroups;
	HDFAtom<string> annotationStringAtom;
	map<string,int> experimentNameToIndex;
	map<int,int>  movieIdToIndex;

	int Initialize(CommonFG &group, string _refGroupName) {
		refGroupName = _refGroupName;
		refGroup.Initialize(group, _refGroupName);
		//		annotationStringAtom.Initialize(refGroup.group, "annotationString");
	}

	HDFCmpExperimentGroup* InitializeExperimentGroup(string experimentGroupName) {
		if (refGroup.ContainsObject(experimentGroupName)) {
			HDFCmpExperimentGroup* newGroup = new HDFCmpExperimentGroup;

			if (newGroup->Initialize(refGroup, experimentGroupName) == 0) {
				cout << "Could not initialize the exp group." << endl;
				exit(0);
			}
			experimentNameToIndex[experimentGroupName] = readGroups.size();
			readGroups.push_back(newGroup);
			return newGroup;
		}
		else {
			return NULL;
		}
	}
};


#endif
