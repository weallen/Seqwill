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

#ifndef DATA_HDF_HDF_BAS_WRITER_H_
#define DATA_HDF_HDF_BAS_WRITER_H_

#include "data/hdf/HDFArray.h"
#include "data/hdf/BufferedHDFArray.h"
#include "data/hdf/HDF2DArray.h"
#include "data/hdf/BufferedHDF2DArray.h"
#include "data/hdf/HDFAtom.h"
#include "data/hdf/HDFFile.h"
#include "DatasetCollection.h"
#include "../../Enumerations.h"
#include "../../utils/SMRTReadUtils.h"
#include "../../FASTQSequence.h"
#include "../../Types.h"
#include <sstream>

using namespace H5;
using namespace std;

class HDFBasWriter : public DatasetCollection {
	HDFFile outFile;
	string hdfFileName;
	string movieName, runCode;
	static const int bufferSize = 16;
	
	float frameRate;
	float numFrames;

	BufferedHDFArray<int> nElemArray;
	BufferedHDFArray<int> zmwXCoordArray;
	BufferedHDFArray<int> zmwYCoordArray;
	BufferedHDFArray<char> baseArray;
	BufferedHDFArray<unsigned char> qualArray;
	BufferedHDFArray<unsigned int> simulatedCoordinateArray;
	BufferedHDFArray<unsigned int> simulatedSequenceIndexArray;

	HDFAtom<string> movieNameAtom, runCodeAtom, platformNameAtom;
	HDFAtom<unsigned int> platformIdAtom;

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
	BufferedHDFArray<HalfWord> preBaseFramesArray;
	BufferedHDFArray<HalfWord> widthInFramesArray;
	BufferedHDFArray<HalfWord> pulseIndexArray;
	BufferedHDF2DArray<unsigned char> preBaseDeletionQVArray;
	HDFGroup rootGroup;
	HDFGroup scanDataGroup;
	HDFGroup runInfoGroup, acqParamsGroup;
	HDFGroup pulseDataGroup;
	HDFGroup baseCallGroup, pulseCallGroup;
	HDFGroup zmwGroup;
	HDFGroup plsZMWGroup;
	PlatformType platformId;
	string   platformName;
 public:
	~HDFBasWriter() { 
		//
		// Assume that flushing out and closing the hdf file must be one
		// manually and not in a destructor.
		//
	}
	void InitializeDefaultIncludedFields() {
		IncludeField("Basecall");
		IncludeField("DeletionQV");
		IncludeField("DeletionTag");
		IncludeField("InsertionQV");
		IncludeField("SubstitutionTag");
		IncludeField("SubstitutionQV");
		IncludeField("QualityValue");
		IncludeField("WidthInFrames");
		IncludeField("PulseIndex");
		IncludeField("PreBaseFrames");
	}
	

	void Flush() {
		nElemArray.Flush();
		if (includedFields["zmwXCoord"])
		zmwXCoordArray.Flush();
		if (includedFields["zmwYCoord"])
			zmwYCoordArray.Flush();
		if (includedFields["Basecall"])
			baseArray.Flush();
		if (includedFields["QualityValue"])
			qualArray.Flush();
		if (includedFields["DeletionQV"])
			deletionQVArray.Flush();
		if (includedFields["DeletionTag"])
			deletionTagArray.Flush();
		if (includedFields["InsertionQV"])
			insertionQVArray.Flush();
		if (includedFields["SubstitutionTag"])
			substitutionTagArray.Flush();
		if (includedFields["SubstitutionQV"])
			substitutionQVArray.Flush();
		if (includedFields["HoleNumber"])
			holeNumberArray.Flush();
		if (includedFields["PreBaseFrames"])
			preBaseFramesArray.Flush();
		if (includedFields["WidthInFrames"])
			widthInFramesArray.Flush();
		if (includedFields["HoleXY"])
			holeXY2D.Flush();
		if (includedFields["SimulatedCoordinate"])
			simulatedCoordinateArray.Flush();
		if (includedFields["SimulatedSequenceIndex"]) 
			simulatedSequenceIndexArray.Flush();
	}

	HDFBasWriter() {
		/*
		 * Default to astro for now.  This may need to change to a NO_ID
		 * platform, in which case it must be set with Initialize().
		 */
		frameRate = 100.0;
		numFrames = 10000000;
		fieldNames.push_back("zmwXCoord");
		fieldNames.push_back("zmwYCoord");
		fieldNames.push_back("QualityValue");
		fieldNames.push_back("Basecall");
		fieldNames.push_back("DeletionQV");
		fieldNames.push_back("DeletionTag");
		fieldNames.push_back("InsertionQV");
		fieldNames.push_back("SubstitutionQV");
		fieldNames.push_back("SubstitutionTag");
		fieldNames.push_back("WidthInFrames");
		fieldNames.push_back("HoleNumber");
		fieldNames.push_back("HoleXY");
		fieldNames.push_back("PreBaseFrames");
		fieldNames.push_back("PulseIndex");
		fieldNames.push_back("SimulatedCoordinate");
		fieldNames.push_back("SimulatedSequenceIndex");
		InitializeAllFields(false);
		platformId = Springfield;
	}

	void Close() {
		Flush();
		outFile.Close();
	}

	void SetPlatform(PlatformType _platform) {
		platformId = _platform;
	}
	void AddMovieName(string movieName) {
		movieNameAtom.Create(runInfoGroup.group, "MovieName",movieName);
	}
	/*
	 * Initialization without a runCode is implicitly a springfield
	 * platform.  You can change it if you really want.
	 */
	void Initialize(string _hdfFileName, string movieName, PlatformType _platform = Springfield) {
		Initialize(_hdfFileName, _platform);
		AddMovieName(movieName);
	}
	void Initialize(string _hdfFileName, string movieName, string runCode, PlatformType _platform = Astro) {
		Initialize(_hdfFileName, _platform);
		if (movieName != "" and runCode != "")
			AddRunInfo(movieName, runCode);
	}

	void AddRunInfo(string movieName, string runCode) {
		AddMovieName(movieName);
		runCodeAtom.Create(runInfoGroup.group, "RunCode", runCode);
	}
	
	void AddPlatformName(string platformName) {
		platformNameAtom.Create(runInfoGroup.group, "PlatformName", platformName);
	}

	void AddPlatformId(PlatformType _platformId) {
		platformIdAtom.Create(runInfoGroup.group, "PlatformId");
		platformIdAtom.Write((unsigned int)_platformId);
	}

	void WriteSimulatedCoordinate(unsigned int coord) {
		simulatedCoordinateArray.Write(&coord,1);
	}

	void WriteSimulatedSequenceIndex(unsigned int index) {
		simulatedSequenceIndexArray.Write(&index,1);
	}

	void Initialize(string _hdfFileName, PlatformType _platform) {
		hdfFileName = _hdfFileName;
		platformId  = _platform;
		outFile.Create(hdfFileName);
		rootGroup.Initialize(*outFile.hdfFile, "/");		
		rootGroup.AddGroup("PulseData"); 
		rootGroup.AddGroup("PulseData/BaseCalls"); 
		rootGroup.AddGroup("PulseData/BaseCalls/ZMW"); 
		rootGroup.AddGroup("ScanData"); 
		rootGroup.AddGroup("ScanData/RunInfo"); 
		rootGroup.AddGroup("ScanData/AcqParams");
		rootGroup.AddGroup("PulseData/PulseCalls");
		rootGroup.AddGroup("PulseData/PulseCalls/ZMW");
		
		scanDataGroup.Initialize(rootGroup.group, "ScanData");
		acqParamsGroup.Initialize(scanDataGroup.group, "AcqParams");


		HDFAtom<float> frameRateAtom;
		HDFAtom<unsigned int> numFramesAtom;
		frameRateAtom.Create(acqParamsGroup.group, "FrameRate");
		numFramesAtom.Create(acqParamsGroup.group, "NumFrames");
		frameRateAtom.Write(frameRate);
		numFramesAtom.Write(numFrames);
		
		pulseDataGroup.Initialize(rootGroup.group, "PulseData");
		baseCallGroup.Initialize(pulseDataGroup.group, "BaseCalls");
		pulseCallGroup.Initialize(pulseDataGroup.group, "PulseCalls");
		runInfoGroup.Initialize(scanDataGroup.group, "RunInfo");
		zmwGroup.Initialize(baseCallGroup.group, "ZMW");

		nElemArray.Initialize(&zmwGroup.group, "NumEvent", bufferSize);
		if (includedFields["Basecall"])
			baseArray.Initialize(&baseCallGroup.group, "Basecall", bufferSize);
		if (includedFields["QualityValue"])
			qualArray.Initialize(&baseCallGroup.group, "QualityValue", bufferSize);
		if (includedFields["DeletionQV"])
			deletionQVArray.Initialize(&baseCallGroup.group, "DeletionQV", bufferSize);
		if (includedFields["DeletionTag"])
			deletionTagArray.Initialize(&baseCallGroup.group, "DeletionTag", bufferSize);
		if (includedFields["InsertionQV"])
			insertionQVArray.Initialize(&baseCallGroup.group, "InsertionQV", bufferSize);
		if (includedFields["PreBaseDeletionQV"])
			preBaseDeletionQVArray.Initialize(&baseCallGroup.group, "PreBaseDeletionQV", 4, bufferSize);
		if (includedFields["SubstitutionTag"])
			substitutionTagArray.Initialize(&baseCallGroup.group,   "SubstituionTag", bufferSize);
		if (includedFields["SubstitutionQV"])
			substitutionQVArray.Initialize(&baseCallGroup.group,     "SubstitutionQV", bufferSize);
		if (includedFields["WidthInFrames"])
			widthInFramesArray.Initialize(&baseCallGroup.group, "WidthInFrames", bufferSize);
		if (includedFields["PreBaseFrames"])
			preBaseFramesArray.Initialize(&baseCallGroup.group, "PreBaseFrames", bufferSize);
		if (includedFields["SimulatedCoordinate"]) {
			simulatedCoordinateArray.Initialize(&baseCallGroup.group, "SimulatedCoordinate", bufferSize);
		}
		if (includedFields["SimulatedSequenceIndex"]) {
			simulatedSequenceIndexArray.Initialize(&baseCallGroup.group, "SimulatedSequenceIndex", bufferSize);
		}

		if (platformId == Astro) {
			holeXY2D.Initialize(&zmwGroup.group, "HoleXY", 2, bufferSize);
		}
		holeNumberArray.Initialize(&zmwGroup.group, "HoleNumber", bufferSize);
	}

	int WriteIdentifiers(UInt holeNumber, int x=0, int y=0 ) {
		//
		// Write hole number regardless of platform type.
		//
		holeNumberArray.Write(&holeNumber, 1);

		if (platformId == Astro) {
			// now extract the x an y coordinates.
			uint16_t xy[2] = {(uint16_t) x, (uint16_t) y};
			holeXY2D.WriteRow(xy, 2);
		}
		return 1;
	}

	int WriteQualities(FASTQSequence &seq) {
		qualArray.Write(seq.qual, seq.length);

		if (includedFields["DeletionQV"] and seq.deletionQV != NULL) {
			deletionQVArray.Write(seq.deletionQV, seq.length);
		}
		if (includedFields["PreBaseDeletionQV"] and seq.preBaseDeletionQV != NULL) {
			DNALength readPos;
			for (readPos = 0; readPos < seq.length; readPos++) {
				preBaseDeletionQVArray.WriteRow(&seq.preBaseDeletionQV[readPos*4], 4);
			}
		}
		if (includedFields["DeletionTag"] and seq.deletionTag != NULL) {
			deletionTagArray.Write(seq.deletionTag, seq.length);
		}
		if (includedFields["InsertionQV"] and seq.insertionQV != NULL) {
			insertionQVArray.Write(seq.insertionQV, seq.length);
		}
		if (includedFields["SubstitutionQV"] and seq.substitutionQV != NULL) {
			substitutionQVArray.Write(seq.substitutionQV, seq.length);
		}
		if (includedFields["SubstitutionTag"] and seq.substitutionTag != NULL) {
			substitutionTagArray.Write(seq.substitutionTag, seq.length);
		}
	}

	int WriteBases(FASTASequence &seq ) {
		int lenArray[1] = {seq.length};
		nElemArray.Write(lenArray, 1);
		baseArray.Write((const char*) seq.seq, seq.length);
		return 1;
	}

	int Write(SMRTSequence &seq) {
		WriteBases(seq);
		WriteQualities(seq);

		if (includedFields["PreBaseFrames"] and seq.preBaseFrames != NULL) {
			preBaseFramesArray.Write(seq.preBaseFrames, seq.length);
		}
		if (includedFields["WidthInFrames"] and seq.widthInFrames != NULL) {
			widthInFramesArray.Write(seq.widthInFrames, seq.length);
		}
		WriteIdentifiers(seq.zmwData.holeNumber, seq.xy[0], seq.xy[1]);
		return 1;
	}

	int Write(FASTQSequence &seq) {

		int x, y;
		UInt holeNumber;

		WriteBases(seq);
		WriteQualities(seq);

		if (platformId == Astro) {
			//
			// now extract the x an y coordinates.
			GetSMRTReadCoordinates(seq, x, y);
			UInt holeNumber;
			seq.GetHoleNumber((int&) holeNumber);
		}
		if( platformId == Springfield){ 
			GetSpringfieldHoleNumberFromTitle(seq, holeNumber);
		}

		WriteIdentifiers(holeNumber,x,y);

		return 1;
	}
};





#endif
