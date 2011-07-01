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

#ifndef DATA_HDF_HDF_GROUP_H_
#define DATA_HDF_HDF_GROUP_H_

#include <H5Cpp.h>
#include "hdf5/HDFAttributable.h"
#include "base/Types.h"
#include "base/StringUtil.h"
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace H5;
using namespace std;

class HDFGroup : public HDFAttributable {
 public:
	vector<string> objectNames;
	string objectName;
	Group group;
	bool  groupIsInitialized;

 HDFGroup() : HDFAttributable() {
		groupIsInitialized = false;
	}

	void AddGroup(string groupName) {
		/*
		 * Add possibly a nested group to the hdf file.
		 */
		group.createGroup(groupName);
		return;

		vector<string> groupPath;
		int nGroups;
        vector<char> sep;
        sep.push_back('/');
		nGroups = Tokenize(groupName, sep, groupPath);
		VectorIndex groupPathIndex;
		Group curGroup;
		string groupSubpathName;
		for (groupPathIndex = 0; groupPathIndex < groupPath.size(); groupPathIndex++) {
      //
      //  Look to see if the curFG contains this part of the path.
      //
			bool groupFound = false;
			groupSubpathName.append(groupPath[groupPathIndex]);
			groupSubpathName.append("/");
			group.createGroup(groupPath[groupPathIndex].c_str());
		}
	}

	int Initialize(CommonFG &fg, string groupName){ 
		try {
			group   = fg.openGroup(groupName.c_str());
			//
			// Get a list of everything in the group, just for kicks.
			//
			hsize_t numGroupObjs;
			numGroupObjs = group.getNumObjs();
			hsize_t objIdx;
			for (objIdx = 0; objIdx < numGroupObjs; objIdx++) {
				H5std_string objName;
				size_t objNameSize;
				objName = group.getObjnameByIdx(objIdx);
				objectNames.push_back(objName);
			}

			StoreAttributeNames(group);
			groupIsInitialized = true;
			return 1;
		}
		catch(FileIException &e) {
			return 0;
		}
		catch(GroupIException &e) {
			return 0;
		}
		catch(Exception &e ) {
			return 0;
		}
		return 1;
	}
	
		
	bool ContainsObject(string objectName) {
		int i;
		for (i = 0; i < objectNames.size(); i++) {
			if (objectNames[i] == objectName) { return true; }
		}
		return false;
	}

	void AddObjectName(string objectName) {
		objectNames.push_back(objectName);
	}

	template<typename T_Object>
	void InitializeObject(T_Object &object, string objectName) {
		if (!ContainsObject(objectName)) {
			object.Initialize(&group, objectName);
			AddObjectName(objectName);
		}
		else {
			//
			// The dataset exists, just open it so it will not be created on
			// a flush later on. 
			//
			object.Initialize(&group, objectName);
			if (object.InitializeDataset(group,objectName) == 0) {
				cout << "ERROR opening " << objectName << endl;
				exit(0);
			}
		}
	}

	template<typename T_Object>
		void Initialize2DObject(T_Object &object, string objectName, int rowLength) {
		object.Initialize(&group, objectName, rowLength);
		AddObjectName(objectName);
	}

	void Close() {
		if (groupIsInitialized) {
			group.close();
		}
	}
};


#endif
