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

#ifndef DATA_HDF_HDF_CMP_DATA_H_
#define DATA_HDF_HDF_CMP_DATA_H_

#include "data/hdf/HDFAtom.h"
#include "data/hdf/HDFCmpRefAlignmentGroup.h"

class HDFCmpData {
 public:
	HDFAtom<string> commandLine;
	HDFAtom<unsigned int> lastAlignmentRow;
	H5File hdfCmpFile;
	HDFArray<int> movieNameIdArray;
	HDFAtom<int>  movieNameIdLastRow;
	HDFArray<string> movieNameArray;

	HDFArray<string> readGroupPathArray;
	HDFArray<int>    readGroupPathIdArray;
	HDFAtom<int>     readGroupPathIdLastRow;

	HDFArray<int>    refSeqNameIdArray;
	HDFAtom<int>     refSeqNameIdLastRow;
	HDFArray<string> refSeqNameArray;
	static const int NCols=22;
	vector<HDFAtom<string> >  colNameAtoms;
	vector<HDFCmpRefAlignmentGroup*> refAlignGroups;
	/*	static const char *colNameIds[] = {"00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
															"10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
															"20", "21"};*/
	static const char *colNameIds[];
	
	void Close() {
		hdfCmpFile.close();
	}
};

const char * HDFCmpData::colNameIds[] = {"00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
																				 "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
																				 "20", "21"};


#endif
