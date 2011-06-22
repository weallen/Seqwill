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

#ifndef DATA_HDF_HDF_CMP_ROOT_GROUP_H_
#define DATA_HDF_HDF_CMP_ROOT_GROUP_H_

#include "HDFAtom.h"
#include "HDF2DArray.h"
#include "datastructures/alignment/CmpFile.h"
template <typename T_Alignment>
class HDFCmpRootGroup {
 public:
	HDFGroup rootGroup;
	HDFAtom<string> version;
	HDFAtom<string> index;
	HDFAtom<string> readType;
	HDFAtom<string> commandLine;
	HDF2DArray<string> fileLog;

	~HDFCmpRootGroup() {
		rootGroup.Close();
	}

	int Initialize(H5File &cmpFile) {
		if (rootGroup.Initialize(cmpFile, "/") == 0) { return 0; }
		if (rootGroup.ContainsObject("Version")) {
			if (version.Initialize(rootGroup.group, "Version") == 0) { return 0; }		
		}
		if (rootGroup.ContainsObject("Index")) {
			if (index.Initialize(rootGroup.group, "Index") == 0) { return 0; }
		}
		if (rootGroup.ContainsObject("ReadType")) {
			if (readType.Initialize(rootGroup.group, "ReadType") == 0) { return 0; }
		}
		if (rootGroup.ContainsObject("CommandLine")) {
			if (commandLine.Initialize(rootGroup.group, "CommandLine") == 0) { return 0; }
		}
		
		//
		// For now, disable file log initialization until
		// hdf2darray<string> is tested more thoroughly 
		//
		// if (fileLog.Initialize(rootGroup.group, "FileLog") == 0) {
		// return 0;}
		//
		return 1;
	}
	
	void ReadAttributes(CmpFile<T_Alignment> &cmpFile) {
		if (rootGroup.ContainsObject("Version")) {
		version.Read(cmpFile.version);
		}
		if (rootGroup.ContainsObject("Index")) {
			index.Read(cmpFile.index);
		}
		if (rootGroup.ContainsObject("ReadType")) {
		readType.Read(cmpFile.readType);
		}
		if (rootGroup.ContainsObject("CommandLine")) {
		commandLine.Read(cmpFile.commandLine);
		}
	}
};
 
	

#endif
