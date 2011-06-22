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

#ifndef DATA_HDF_HDF_PLS_READER_H_
#define DATA_HDF_HDF_PLS_READER_H_

#include "data/hdf/HDFArray.h"
#include "data/hdf/HDF2DArray.h"
#include "data/hdf/HDFAtom.h"
#include "data/hdf/PlatformId.h"
#include "data/hdf/HDFGroup.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFPulseDataFile.h"
#include "data/hdf/DatasetCollection.h"
#include "data/hdf/HDFZMWReader.h"
#include "datastructures/reads/PulseFile.h"
#include "HDFScanDataReader.h"
#include "FASTQSequence.h"
#include <sstream>
#include <vector>
#include <set>
using namespace H5;
using namespace std;
/*
 * Interface for reading pulse information from a .pls.h5 file.
 * To read both pls and bas information, use the HDFBasReader.
 */

class HDFPlsReader : public DatasetCollection, public HDFPulseDataFile  {
	DNALength curBasePos;
	int curRead;
	int nReads;

	HDFGroup pulseCallsGroup;

	HDF2DArray<uint16_t> meanSignalArray;
	HDF2DArray<uint16_t> midSignalArray;
	HDF2DArray<uint16_t> maxSignalArray;
	HDFArray<uint16_t>   plsWidthInFramesArray;
	HDFArray<unsigned int> startFrameArray;
	HDFArray<float>      classifierQVArray;

	HDFZMWReader zmwReader;

 public:
 HDFPlsReader() : HDFPulseDataFile(), DatasetCollection() {
		fieldNames.push_back("MeanSignal");
		fieldNames.push_back("MidSignal");
		fieldNames.push_back("MaxSignal");
		fieldNames.push_back("StartFrame");
		fieldNames.push_back("ClassifierQV");
		fieldNames.push_back("WidthInFrames");
		fieldNames.push_back("FrameRate");
		fieldNames.push_back("WhenStarted");
		fieldNames.push_back("NumEvent"); // this is read from the zmw, but control it here.
		InitializeAllFields(false);
	}

	int InitializeCommon() {
		return 1;
	}

	int Initialize(string hdfPlsFileName) {
		/*
		 * Initialize access to the HDF file.
		 */

		if (HDFPulseDataFile::Initialize(hdfPlsFileName) == 0) {
			return 0;
		}

		if (pulseDataGroup.ContainsObject("PulseCalls") == 0 or
				pulseCallsGroup.Initialize(pulseDataGroup.group, "PulseCalls") == 0) {
			return 0;
		}
		zmwReader.Initialize(&pulseCallsGroup);


		//
		// Initialize arrays for reading.  This has been preconfigured to 
		//
		if (!InitializePulseGroup()) return 0;
		if (!InitializeDataset(pulseCallsGroup, meanSignalArray,      "MeanSignal"))    return 0;
		if (!InitializeDataset(pulseCallsGroup, midSignalArray,       "MidSignal"))     return 0;
		if (!InitializeDataset(pulseCallsGroup, maxSignalArray,       "MaxSignal"))     return 0;
		if (!InitializeDataset(pulseCallsGroup, startFrameArray,      "StartFrame"))    return 0;
		if (!InitializeDataset(pulseCallsGroup, plsWidthInFramesArray,"WidthInFrames")) return 0;
		if (!InitializeDataset(pulseCallsGroup, classifierQVArray,    "ClassifierQV"))  return 0;
		return 1;
	}

	void GetAllMeanSignal(vector<uint16_t> &meanSignal) {
		CheckMemoryAllocation(meanSignalArray.GetNCols() * meanSignalArray.GetNRows(), maxAllocNElements, "MeanSignal");
		meanSignal.resize(meanSignalArray.GetNCols() * meanSignalArray.GetNRows());
		meanSignalArray.Read(0, meanSignalArray.GetNRows(), &meanSignal[0]);
	}

	void GetAllMidSignal(vector<uint16_t> &midSignal) {
		CheckMemoryAllocation(midSignalArray.GetNCols() * midSignalArray.GetNRows(), maxAllocNElements, "MidSignal");
		midSignal.resize(midSignalArray.GetNCols() * midSignalArray.GetNRows());
		midSignalArray.Read(0, midSignalArray.GetNRows(), &midSignal[0]);
	}

	void GetAllMaxSignal(vector<uint16_t> &maxSignal) {
		CheckMemoryAllocation(maxSignalArray.GetNCols() * maxSignalArray.GetNRows(), maxAllocNElements, "MaxSignal");
		maxSignal.resize(maxSignalArray.GetNCols() * maxSignalArray.GetNRows());
		maxSignalArray.Read(0, maxSignalArray.GetNRows(), &maxSignal[0]);
	}

	void GetAllStartFrames(vector<UInt> &startFrame) {
		CheckMemoryAllocation(startFrameArray.arrayLength, maxAllocNElements, "StartFrame");
		startFrame.resize(startFrameArray.arrayLength);
		startFrameArray.Read(0,startFrameArray.arrayLength, &startFrame[0]);
	}

	void GetAllPlsWidthInFrames(vector<uint16_t> &widthInFrames) {
		CheckMemoryAllocation(plsWidthInFramesArray.arrayLength, maxAllocNElements, "WidthInFrames (pulse)");
		widthInFrames.resize(plsWidthInFramesArray.arrayLength);
		plsWidthInFramesArray.Read(0,plsWidthInFramesArray.arrayLength, &widthInFrames[0]);
	}
	
	void GetAllClassifierQV(vector<float> &classifierQV) {
		CheckMemoryAllocation(classifierQVArray.arrayLength, maxAllocNElements, "ClassifierQV (pulse)");
		classifierQV.resize(classifierQVArray.arrayLength);
		classifierQVArray.Read(0, classifierQVArray.arrayLength, &classifierQV[0]);
	}

	void GetAllNumEvent(vector<int> &numEvent) {
		CheckMemoryAllocation(zmwReader.numEventArray.arrayLength, maxAllocNElements, "NumEvent (pulse)");
		numEvent.resize(zmwReader.numEventArray.arrayLength);
		zmwReader.numEventArray.Read(0, zmwReader.numEventArray.arrayLength, &numEvent[0]);
	}

	void ReadPulseFile(PulseFile &pulseFile) {
		if (scanDataReader.fileHasScanData) {
			scanDataReader.Read(pulseFile.scanData);
		}
		if (includedFields["StartFrame"]) {
			GetAllStartFrames(pulseFile.startFrame);
		}
		if (includedFields["WidthInFrames"]) {
			GetAllPlsWidthInFrames(pulseFile.plsWidthInFrames);
		}
		if (includedFields["MeanSignal"]) {
			GetAllMeanSignal(pulseFile.meanSignal);
		}
		if (includedFields["MidSignal"]) {
			GetAllMidSignal(pulseFile.midSignal);
		}
		if (includedFields["MaxSignal"]) {
			GetAllMaxSignal(pulseFile.maxSignal);
		}
		if (includedFields["NumEvent"]) {
			GetAllNumEvent(pulseFile.numEvent);
		}
		if (includedFields["ClassifierQV"]) {
			GetAllClassifierQV(pulseFile.classifierQV);
		}
	}

	int GetNext() {
		return 1;
	}
	

	void Close() {
		
	}
};


#endif
