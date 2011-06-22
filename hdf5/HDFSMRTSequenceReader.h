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

#ifndef DATA_HDF_HDF_SMRT_SEQUENCE_READER
#define DATA_HDF_HDF_SMRT_SEQUENCE_READER

#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFZMWReader.h"
#include "SMRTSequence.h"

template<typename T_SMRT_Sequence>
class HDFSMRTSequenceReader : public HDFBasReader {
 public:
	HDFZMWReader zmwReader;
	bool readQuality;
	int Initialize(string hdfBasFileName, bool _readQuality=true) {
		HDFBasReader::Initialize(hdfBasFileName);
		zmwReader.Initialize(hdfBasFile);
		readQuality = _readQuality;
		if (baseCallsGroup.ContainsObject("WidthInFrames")) {
			useWidthInFrames = InitializeField(baseCallsGroup, "WidthInFrames");
		}
		if (baseCallsGroup.ContainsObject("PreBaseFrames")) {
			usePreBaseFrames = InitializeField(baseCallsGroup, "PreBaseFrames");
		}
		if (baseCallsGroup.ContainsObject("PulseIndex")) {
			usePulseIndex = InitializeField(baseCallsGroup, "PulseIndex");
		}
		
	}


	int GetNext(T_SMRT_Sequence &seq) {
		int retVal;
		//
		// Copy the curBasePos from the bas reader since it gets advanced
		// in the GetNext function.
		//
		DNALength curBasePosCopy = curBasePos;
		retVal = HDFBasReader::GetNext(seq);
		//
		// Bail now if the file is already done
		if (retVal == 0) {
			return retVal;
		}
		zmwReader.GetNext(seq.zmwData);
		return retVal;
	}
	int Advance(int nSteps) {
		int retVal;
		retVal = HDFBasReader::Advance(nSteps);
		zmwReader.Advance(nSteps);
		return retVal;
	}


};

template<>
int HDFSMRTSequenceReader<FASTASequence>::GetNext(FASTASequence &seq) {
		int retVal;
		if (readQuality) {
			retVal = HDFBasReader::GetNext(seq);
		}
		//
		// Bail now if the file is already done
		if (retVal == 0) {
			return retVal;
		}
		zmwReader.GetNext(seq.zmwData);
		return retVal;
	}

#endif
