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

#ifndef DATA_HDF_HDF_PLS_WRITER_H_
#define DATA_HDF_HDF_PLS_WRITER_H_

#include "data/hdf/HDFArray.h"
#include "data/hdf/BufferedHDFArray.h"
#include "data/hdf/HDF2DArray.h"
#include "data/hdf/BufferedHDF2DArray.h"
#include "data/hdf/HDFAtom.h"
#include "data/hdf/HDFFile.h"
#include "data/hdf/PlatformId.h"
#include "utils/SMRTReadUtils.h"
#include "FASTQSequence.h"
#include <sstream>

using namespace H5;
using namespace std;

class HDFPlsWriter {
	HDFFile outFile;
	string hdfFileName;
	string movieName, runCode;
	PlatformId platformId;
	static const int bufferSize = 16;
	
	BufferedHDFArray<int> nElemArray;
	BufferedHDFArray<int> zmwXCoordArray;
	BufferedHDFArray<int> zmwYCoordArray;
	BufferedHDFArray<char> baseArray;
	BufferedHDFArray<unsigned char> qualArray;

	/*
	HDFArray<int> nElemArray;
	HDFArray<int> zmwXCoordArray;
	HDFArray<int> zmwYCoordArray;
	HDFArray<char> baseArray;
	HDFArray<unsigned char> qualArray;
	*/
	HDFAtom<string> movieNameAtom, runCodeAtom;

	//
	// Astro specific arrays.
	//
	BufferedHDF2DArray<uint16_t> holeXY2D;

	//
	// Springfield specific arrays.
	//
	BufferedHDFArray<unsigned int> holeNumberArray;
	
	//
	// Define arrays for rich quality values.
	// 

	BufferedHDFArray<unsigned char> deletionQVArray;
	BufferedHDFArray<unsigned char> deletionTagArray;
	BufferedHDFArray<unsigned char> insertionQVArray;
	BufferedHDFArray<unsigned char> substitutionTagArray;
	BufferedHDFArray<unsigned char> substitutionQVArray;


	BufferedHDF2DArray<unsigned char> preBaseDeletionQVArray;
	HDFGroup rootGroup;
	Group runInfoGroup;
	Group baseCallGroup;
	Group zmwGroup;
 public:
	~HDFPlsWriter() { 
		nElemArray.Flush();
		zmwXCoordArray.Flush();
		zmwYCoordArray.Flush();
		baseArray.Flush();
		qualArray.Flush();
		deletionQVArray.Flush();
		deletionTagArray.Flush();
		insertionQVArray.Flush();
		substitutionTagArray.Flush();
		substitutionQVArray.Flush();
		holeNumberArray.Flush();
	}
	
	HDFPlsWriter() {
		/*
		 * Default to astro for now.  This may need to change to a NO_ID
		 * platform, in which case it must be set with Initialize().
		 */
		platformId = AstroPlatform;
	}
	void AddMovieName(string movieName) {
		movieNameAtom.Create(runInfoGroup, "MovieName",movieName);
	}
	/*
	 * Initialization without a runCode is implicitly a springfield
	 * platform.  You can change it if you really want.
	 */
	void Initialize(string _hdfFileName, string movieName, PlatformId _platformId = SpringfieldPlatform) {
		Initialize(_hdfFileName, _platformId);
		AddMovieName(movieName);
	}
	void Initialize(string _hdfFileName, string movieName, string runCode, PlatformId _platformId = AstroPlatform) {
		Initialize(_hdfFileName, _platformId);
		if (movieName != "" and runCode != "")
			AddRunInfo(movieName, runCode);
	}

	void AddRunInfo(string movieName, string runCode) {
		AddMovieName(movieName);
		runCodeAtom.Create(runInfoGroup, "RunCode", runCode);
	}

	void Initialize(string _hdfFileName, PlatformId _platformId) {
		hdfFileName = _hdfFileName;
		platformId  = _platformId;
		outFile.Create(hdfFileName);
		rootGroup.Initialize(*outFile.hdfFile, "/");
		rootGroup.AddGruop("PulseData");
		rootGroup.AddGroup("PulseData/BaesCalls");
		rootGroup.AddGroup("PulseData/BaseCalls/ZMW"); 
		rootGroup.AddGroup("ScanData/RunInfo"); 
		outFile.OpenGroup("ScanData/RunInfo", runInfoGroup);
		outFile.OpenGroup("PulseData/BaseCalls", baseCallGroup);
		outFile.OpenGroup("PulseData/BaseCalls/ZMW", zmwGroup);

		nElemArray.Initialize(&zmwGroup, "NumEvent", bufferSize);
		baseArray.Initialize(&baseCallGroup, "Basecall", bufferSize);
		qualArray.Initialize(&baseCallGroup, "QualityValue", bufferSize);

		deletionQVArray.Initialize(&baseCallGroup, "DeletionQV", bufferSize);
		deletionTagArray.Initialize(&baseCallGroup, "DeletionTag", bufferSize);
		insertionQVArray.Initialize(&baseCallGroup, "InsertionQV", bufferSize);
		preBaseDeletionQVArray.Initialize(&baseCallGroup, "PreBaseDeletionQV", 4, bufferSize);
		substitutionTagArray.Initialize(&baseCallGroup, "SubstituionTag", bufferSize);
		substitutionQVArray.Initialize(&baseCallGroup, "SubstitutionQV", bufferSize);

		if (platformId == AstroPlatform) {
			holeXY2D.Initialize(&zmwGroup, "HoleXY", 2, bufferSize);
		}
		else if (platformId == SpringfieldPlatform) {
			holeNumberArray.Initialize(&zmwGroup, "HoleNumber", bufferSize);
		}
	}

	int Write(FASTQSequence &seq) {
		int lenArray[1] = {seq.length};
		nElemArray.Write(lenArray, 1);
		qualArray.Write(seq.qual, seq.length);
		baseArray.Write((const char*) seq.seq, seq.length);

		if (seq.deletionQV != NULL) {
			deletionQVArray.Write(seq.deletionQV, seq.length);
		}
		if (seq.preBaseDeletionQV != NULL) {
			DNALength readPos;
			for (readPos = 0; readPos < seq.length; readPos++) {
				preBaseDeletionQVArray.WriteRow(&seq.preBaseDeletionQV[readPos*4], 4);
			}
		}
		if (seq.deletionTag != NULL) {
			deletionTagArray.Write(seq.deletionTag, seq.length);
		}
		if (seq.insertionQV != NULL) {
			insertionQVArray.Write(seq.insertionQV, seq.length);
		}
		if (seq.substitutionQV != NULL) {
			substitutionQVArray.Write(seq.substitutionQV, seq.length);
		}
		if (seq.substitutionTag != NULL) {
			substitutionTagArray.Write(seq.substitutionTag, seq.length);
		}

		if (platformId == AstroPlatform) {
			// now extract the x an y coordinates.
			int x, y;
			GetSMRTReadCoordinates(seq, x, y);
			uint16_t xy[2] = {(uint16_t) x, (uint16_t) y};
			holeXY2D.WriteRow(xy, 2);
			int holeNumber = 0;
			seq.GetHoleNumber(holeNumber);
			holeNumberArray.Write(&holeNumber, 1);			
		}
		else if( platformId == SpringfieldPlatform){ 
			unsigned int holeNumber;
			GetSpringfieldHoleNumberFromTitle(seq, holeNumber);
			holeNumberArray.Write(&holeNumber, 1);
		}
		// For now say this always works. HDF will choke if a problem
		// happens.
		return 1;
	}
};





#endif
