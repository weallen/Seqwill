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

#ifndef HDF_MOVIE_INFO_GROUP_H_
#define HDF_MOVIE_INFO_GROUP_H_

#include "HDFGroup.h"
#include "HDFArray.h"

class HDFMovieInfoGroup {
 public:
	
	HDFGroup movieInfoGroup;
	HDFAtom<string> masterDatasetAtom;
	HDFAtom<UInt>  nRow;
	HDFArray<UInt> idArray;
	HDFArray<UInt> experimentArray;
	HDFArray<string>  nameArray;
	
	HDFArray<string> whenStartedArray;
	
	~HDFMovieInfoGroup() {
		movieInfoGroup.Close();
	}
	int Initialize(HDFGroup &parentGroup) {
		if (movieInfoGroup.Initialize(parentGroup.group, "MovieInfo") == 0) { return 0; }

		if (masterDatasetAtom.Initialize(movieInfoGroup.group, "MasterDataset") == 0) { return 0;}

		if (idArray.Initialize(movieInfoGroup, "ID") == 0) { return 0; }
		if (nameArray.Initialize(movieInfoGroup, "Name") == 0) { return 0; }
		if (movieInfoGroup.ContainsAttribute("nRow")) {
			nRow.Initialize(movieInfoGroup.group, "nRow");
		}
		return 1;
	}
	
	void Read(MovieInfo &movieInfo) {
		int nId = idArray.arrayLength;
		movieInfo.id.resize(nId);
		idArray.Read(0, nId, &movieInfo.id[0]);
		
		int nName = nameArray.arrayLength;
		movieInfo.name.resize(nName);
		int i;
		for (i = 0; i < nName; i++ ){
			nameArray.Read(i,i+1,&movieInfo.name[i]);
		}
	}
};
			
			
			
				
#endif
