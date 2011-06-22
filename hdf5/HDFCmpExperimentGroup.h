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

#ifndef DATA_HDF_HDF_CMP_EXPERIMENT_GROUP_H_
#define DATA_HDF_HDF_CMP_EXPERIMENT_GROUP_H_

#include "data/hdf/HDFGroup.h"
#include "data/hdf/HDFArray.h"
#include "data/hdf/HDFAtom.h"
#include "Types.h"

class HDFCmpExperimentGroup {
 public:
	int Initialize(HDFGroup &refGroup, string experimentGroupName) {
		if (experimentGroup.Initialize(refGroup.group, experimentGroupName) == 0) { return 0; }
		
		if (masterDatasetAtom.Initialize(experimentGroup.group, "MasterDataset") == 0) { return 0; }
		if (experimentGroup.ContainsAttribute("nrow")) {
			nRowAtom.Initialize(experimentGroup.group, "nRow");
		}
		if (alignmentArray.Initialize(experimentGroup, "AlnArray") == 0) { return 0; }

		return 1;
	}

	void Read();
	HDFAtom<string> masterDatasetAtom;
	HDFAtom<UInt>   nRowAtom;
	HDFGroup experimentGroup;
	HDFArray<unsigned char> alignmentArray;
	HDFArray<float> startTimeOffset;
	HDFArray<float> qualityValue;
	HDFArray<float> ipd;
	HDFArray<float> preBaseFrames;
	HDFArray<float> deletionQV;
	HDFArray<float> insertionQV;
	HDFArray<float> classifierQV;
	HDFArray<float> substitutionQV;
	HDFArray<float> light;
	HDFArray<float> widthInFrames;
	HDFArray<float> pulseWidth;
	HDFArray<float> startTime;
	HDFArray<float> pkmid;
	HDFArray<float> pkmax;
	HDFArray<float> pkmean;
	HDFAtom<int>    alignmentArrayLastRow;
};
#endif
