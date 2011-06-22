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

#ifndef DATA_HDF_HDF_BAS_READER_H_
#define DATA_HDF_HDF_BAS_READER_H_

#include "data/hdf/DatasetCollection.h"
#include "data/hdf/HDFArray.h"
#include "data/hdf/HDF2DArray.h"
#include "data/hdf/HDFAtom.h"
#include "data/hdf/HDFGroup.h"
#include "datastructures/reads/BaseFile.h"
#include "FASTQSequence.h"
#include "SMRTSequence.h"
#include "Enumerations.h"
#include "data/hdf/HDFZMWReader.h"
#include "HDFScanDataReader.h"
#include "HDFPulseDataFile.h"
#include <stdlib.h>
#include <sstream>
#include <vector>

using namespace H5;
using namespace std;

/*
 * Below is sample code for using the bas reader to read in sequences
 * from a .bas.h5 or .pls.h5 file.
 *
#include "data/hdf/HDFBasReader.h"
#include "FASTQSequence.h"
#include <string>

int main(int argc, char* argv[]) {
	
	string basFile = argv[1];

	HDFBasReader reader;

	reader.Initialize(basFile);

	FASTQSequence read;

	while(reader.GetNextQ(read)) {
		read.PrintQualSeq(cout);
	}


	return 0;
}

 * If you wanted to read only fasta sequences and not fastq, use:
 
 FASTASequence read;
 while (reader.GetNext(read) {
   read.PritnSeq(cout);
 }

 */

template<typename T_Sequence>
class T_HDFBasReader : public DatasetCollection, public HDFPulseDataFile {
 public:
	DNALength curBasePos;
	int curRead;
	int nReads;
	unsigned int nBases;


	bool readPulseInformation;
	bool hasRegionTable;

	HDFArray<int> nElemArray;
	HDFArray<int> zmwXCoordArray;
	HDFArray<int> zmwYCoordArray;
	HDFArray<unsigned char> baseArray;
	HDFArray<unsigned char> deletionQVArray;
	HDFArray<unsigned char> deletionTagArray;
	HDFArray<unsigned char> insertionQVArray;
	HDFArray<unsigned char> substitutionTagArray;
	HDFArray<unsigned char> substitutionQVArray;
	HDFArray<unsigned char> qualArray;
	HDFArray<unsigned int> simulatedCoordinateArray;
	HDFArray<unsigned int> simulatedSequenceIndexArray;
	HDFArray<uint16_t> basWidthInFramesArray;
	HDFArray<uint16_t> preBaseFramesArray;
	HDFArray<int> pulseIndexArray;
	HDFArray<int> holeIndexArray;
	HDFZMWReader zmwReader;
	HDFGroup baseCallsGroup;
	HDFGroup zmwGroup;
	HDFArray<HalfWord> pulseWidthArray;
	
	bool useWidthInFrames, usePulseWidth, usePulseIndex, useZmwReader;
	
	string baseCallsGroupName;
	bool qualityFieldsAreCritical;
	
	bool useHoleNumbers;
	bool usePreBaseFrames;
	bool useBasHoleXY;
	bool useBasecall;
	bool useQuality;
	
 public:
	PlatformId GetPlatform() {
		return scanDataReader.platformId;
	}

	string GetMovieName() {
		if (scanDataReader.useMovieName) {
			return scanDataReader.movieName;
		}
		else {
			return "";
		}
	}

	string GetRunCode() {
		return scanDataReader.GetRunCode();
	}

	T_HDFBasReader() {
		nReads       = 0;
		curRead      = 0;
		curBasePos   = 0;
		baseCallsGroupName = "BaseCalls";
		qualityFieldsAreCritical = true;
		useZmwReader = false;
		fieldNames.push_back("Basecall");
		fieldNames.push_back("DeletionQV");
		fieldNames.push_back("DeletionTag");
		fieldNames.push_back("InsertionQV");
		fieldNames.push_back("SubstitutionTag");
		fieldNames.push_back("SubstitutionQV");
		fieldNames.push_back("QualityValue");
		fieldNames.push_back("WidthInFrames");
		fieldNames.push_back("PulseIndex");
		fieldNames.push_back("PreBaseFrames");
		fieldNames.push_back("SimulatedCoordinate");
		fieldNames.push_back("SimulatedSequenceIndex");
		InitializeAllFields(false);
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
	
	void InitializeDefaultRequiredFields() {
		requiredFields["Basecall"] = true;
	}

	bool HasRegionTable() {
		return hasRegionTable;
	}
	
	
	bool HasSimulatedCoordinatesStored() {
		if (baseCallsGroup.ContainsObject("SimulatedCoordinate") and
				baseCallsGroup.ContainsObject("SimulatedSequenceIndex")) {
			return true;
		}
		else {
			return false;
		}
	}

	int InitializeForReadingBases() {

		//
		// Initialize root group + scan data information.
		//
		if (HDFPulseDataFile::Initialize(rootGroupPtr) == 0) return 0;

		//
		// Open the base group, this contains all the required information.  
		//
		if (pulseDataGroup.ContainsObject(baseCallsGroupName) == 0 or
				baseCallsGroup.Initialize(pulseDataGroup.group, baseCallsGroupName) == 0) {
			return 0;
		}

		if (pulseDataGroup.ContainsObject("Regions")) {
			hasRegionTable = true;
		}
		else {
			hasRegionTable = false;
		}

		//
		// Initialize read and quality arrays for reading.
		//
		this->InitializeSequenceFields(baseCallsGroup);

		//
		// Initialize simulated coordinate fields if they exist.  They are
		// automatically opened and read from when they are in the bas.h5
		// file since they are used for debugging.
		//
		
		if (baseCallsGroup.ContainsObject("SimulatedCoordinate")) {
			includedFields["SimulatedCoordinate"] = true;
			InitializeDataset(baseCallsGroup, simulatedCoordinateArray, "SimulatedCoordinate");
		}
		else {
			includedFields["SimulatedCoordinate"] = false;
		}

		if (baseCallsGroup.ContainsObject("SimulatedSequenceIndex")) {
			includedFields["SimulatedSequenceIndex"] = true;
			InitializeDataset(baseCallsGroup, simulatedSequenceIndexArray, "SimulatedSequenceIndex");
		}
		else {
			includedFields["SimulatedSequenceIndex"] = false;
		}
		nBases = baseArray.arrayLength;

		return 1;
	}


	int InitializeCommon() {

		//
		// Initialize the smallest set of fields required to import bases.
		//
		if (InitializeForReadingBases() == 0) {
			return 0;
		}
		

		return 1;
	}

	template<typename T_Dataset>
	void InitializeRequiredField(HDFGroup &group, string arrayName, T_Dataset &field) {
		if (group.ContainsObject(arrayName)) {
			if (field.Initialize(group, arrayName) != 0) {
				return;
			}
		}
		cout << "ERROR. Could not initialize dataset " << arrayName << endl;
		exit(0);
	}
		
	template<typename T>
	int InitializeField(HDFGroup &rootGroup, string arrayName,
											T &field, bool &initialized) {
		initialized = false;
		if (rootGroup.ContainsObject(arrayName)) {
			if (field.Initialize(rootGroup, arrayName) != 0) {
				initialized = true;
				return true;
			}
		}
		return false;
	}

	template<typename T>
	int InitializeAttribute(HDFGroup &rootGroup, string arrayName,
											T &field, bool &initialized,
											bool fieldIsCritical = true) {
		int success = 1;
		initialized = false;
		if (rootGroup.ContainsAttribute(arrayName)) {
			if (field.Initialize(rootGroup, arrayName) == 0) {
				success = 0;
			}
			else {
				initialized = true;
			}
		}
		else {
			// the field does not exist
			success = 0;
		}
		if (fieldIsCritical) {
			return success;
		}
		else {
			return 1;
		}
	}
														 
	int InitializeSequenceFields(HDFGroup &baseCallsGroup) {
		//
		// The only field that is absoultely required is 
		if (InitializeDataset(baseCallsGroup, baseArray,             "Basecall")        == false) return 0;
		if (InitializeDataset(baseCallsGroup, qualArray,             "QualityValue")    == false) return 0;
		if (InitializeDataset(baseCallsGroup, insertionQVArray,      "InsertionQV")     == false) return 0;
		if (InitializeDataset(baseCallsGroup, deletionQVArray,       "DeletionQV")      == false) return 0;
		if (InitializeDataset(baseCallsGroup, deletionTagArray,      "DeletionTag")     == false) return 0;
		if (InitializeDataset(baseCallsGroup, substitutionQVArray,   "SubstitutionQV")  == false) return 0;
		if (InitializeDataset(baseCallsGroup, substitutionTagArray,  "SubstitutionTag") == false) return 0;
		if (InitializeDataset(baseCallsGroup, preBaseFramesArray,    "PreBaseFrames")   == false) return 0;
		if (InitializeDataset(baseCallsGroup, pulseIndexArray,       "PulseIndex")      == false) return 0;
		if (InitializeDataset(baseCallsGroup, basWidthInFramesArray, "WidthInFrames")   == false) return 0;

	}

  int InitializeAstro() {
		useBasHoleXY = true;
		return 1;
	}

	int InitializeSpringfield() {
		//
		// For now, no special initialization is required.
		//
		return 1;
	}

	int Initialize(HDFGroup *rootGroupP) {
		rootGroupPtr= rootGroupP;
		return Initialize();
	}

	int Initialize() {
		
		// Return 0 if any of the array inializations do not work.
		if (InitializeCommon() == 0) {
			return 0;
		}

		// Must have zmw information.
		if (zmwReader.Initialize(&baseCallsGroup) == 0) {
			return 0;
		}
		else {
			useZmwReader = true;
		}
		//
		// Get information about the chip - how many zmw's, the hole
		// number indices to keep track of reads, etc..
		//
		nReads = zmwReader.numEventArray.arrayLength;
		


		if (scanDataReader.platformId == AstroPlatform) {
			if (InitializeAstro() == 0) {
				return 0;
			}
		}
		else if (scanDataReader.platformId == SpringfieldPlatform) {
			if (InitializeSpringfield() == 0) {
				return 0;
			}
		}

		/*
		 * Initialize state variables.
		 */
		
		curBasePos = 0;
		curRead    = 0;

		/*
		 * All ok, return that.
		 */
		return 1;
	}

	int Initialize(string hdfBasFileName) {
		/*
		 * Initialize access to the HDF file.  For reading bas files, this
		 * involves:
		 *   - Opening the file, and initializing both base and zmw grops.
		 */
		if (OpenHDFFile(hdfBasFileName) == 0) {
			return 0;
		}

		if (rootGroup.Initialize(hdfBasFile, "/") == 0) {
			return 0;
		}
		rootGroupPtr = &rootGroup;
		return Initialize();
	}

	int GetNumReads() {
		return nReads;
	}

	void AstroBuildReadTitle(string run, string movieTitle, int holeX, int holeY, string &readTitle, unsigned int holeNumber = 0, unsigned int simIndex=0, unsigned int simCoordinate=0) {
		stringstream readTitleStrm;
		if (includedFields["SimulatedSequenceIndex"] == true and
				includedFields["SimulatedCoordinate"] == true) {
			readTitleStrm << "x" << holeX << "_y" << holeY << "_" 
										<< run << "_" << movieTitle
										<< "/holeNumber_" << holeNumber
										<< "/chrIndex_" << simIndex 
										<< "/position_"<< simCoordinate;
		}
		else {
			readTitleStrm << "x" << holeX << "_y" << holeY << "_" << run << "_" << movieTitle;
		}
		readTitle = readTitleStrm.str();
	}

	void SpringfieldBuildReadTitle(string movieTitle, unsigned int holeNumber, string &readTitle, unsigned int simIndex=0, unsigned int simCoordinate=0) {
		stringstream readTitleStrm;
		if (includedFields["SimulatedSequenceIndex"] == true and
				includedFields["SimulatedCoordinate"] == true) {
			readTitleStrm << movieTitle << "/" << holeNumber << "/chrIndex_" << simIndex << "/position_"<< simCoordinate;
		}
		else {
			readTitleStrm << movieTitle << "/" << holeNumber;
		}
		readTitle = readTitleStrm.str();
	}

	int GetNext(FASTASequence &seq) {
		if (curRead == nReads) {
			return 0;
		}
	
		int seqLength;
		seqLength = GetNextWithoutPosAdvance(seq);
		curBasePos += seqLength;
		seq.StorePlatformType(scanDataReader.platformId);
		return 1;
	}

	
	int GetNext(FASTQSequence &seq) {
		if (curRead == nReads) {
			return 0;
		}
		int seqLength = GetNextWithoutPosAdvance(seq);
		seq.length = seqLength;
		if (seqLength > 0 ) {
			if (includedFields["QualityValue"]) {
				seq.AllocateQualitySpace(seqLength);
				qualArray.Read((int)curBasePos, (int) curBasePos + seqLength, (unsigned char*) seq.qual);
			}
		}

		if (includedFields["DeletionQV"]) {
			GetNextDeletionQV(seq);
		}
		if (includedFields["DeletionTag"]) {
			GetNextDeletionTag(seq);
		}
		if (includedFields["InsertionQV"]) {
			GetNextInsertionQV(seq);
		}
		if (includedFields["SubstitutionQV"]) {
			GetNextSubstitutionQV(seq);
		}
		if (includedFields["SubstitutionTag"]) {
			GetNextSubstitutionTag(seq);
		}


		curBasePos += seqLength;
		return 1;
	}

//
// Reading of SMRT Sequences reads both the sequence fields, and the
// fields with ZMW information for identification of this read.
//

 int GetNext(SMRTSequence &seq) {
	 //
	 // Read in quality values.
	 //
	 int retVal;
	 
	 DNALength  curBasPosCopy = curBasePos;
	 //
	 // Getting next advances the curBasPos to the end of 
	 // the current sequence. 
	 //
	 retVal = this->GetNext((FASTQSequence&)seq);
	 if (retVal  == 0) {
		 return 0;
	 }

	 DNALength nextBasePos = curBasePos;
	 curBasePos = curBasPosCopy;

	 if (includedFields["WidthInFrames"]) {
		 assert(nextBasePos <= basWidthInFramesArray.arrayLength);
		 GetNextWidthInFrames(seq);
	 }
	 if (includedFields["PreBaseFrames"]) { 
		 GetNextPreBaseFrames(seq);
	 }
	 if (includedFields["PulseIndex"]) { 
		 GetNextPulseIndex(seq);
	 }
	 curBasePos = nextBasePos;
	 
	 //
	 // Bail now if the file is already done
	 if (retVal == 0) {
		 return retVal;
	 }
	 //
	 // By default, the subread of a read without subread information is
	 // the whole read.
	 //
	 seq.subreadStart = 0;
	 seq.subreadEnd   = seq.length;
	 zmwReader.GetNext(seq.zmwData);
	 return retVal;
 }
 
	int GetAllReadLengths(vector<DNALength> &readLengths) {
		readLengths.resize(nReads);
		nElemArray.Read(0,nReads, (int*) &readLengths[0]);
		return readLengths.size();
	}

	void GetAllPulseIndex(vector<int> &pulseIndex) {
		CheckMemoryAllocation(pulseIndexArray.arrayLength, maxAllocNElements, "PulseIndex");
		pulseIndex.resize(pulseIndexArray.arrayLength);
		pulseIndexArray.Read(0, pulseIndexArray.arrayLength, &pulseIndex[0]);
	}

	int GetAllPreBaseFrames(vector<uint16_t> &preBaseFrames) {
		CheckMemoryAllocation(preBaseFramesArray.arrayLength, maxAllocNElements, "PreBaseFrames");
		preBaseFrames.resize(nBases);
		preBaseFramesArray.Read(0, nBases, &preBaseFrames[0]);
	}

	int GetAllWidthInFrames(vector<uint16_t> &widthInFrames) { 
		CheckMemoryAllocation(basWidthInFramesArray.arrayLength, maxAllocNElements, "WidthInFrames");
		widthInFrames.resize(nBases);
		basWidthInFramesArray.Read(0, nBases, &widthInFrames[0]);
	}

	int GetAllHoleNumbers(vector<unsigned int> &holeNumbers) {
		CheckMemoryAllocation(zmwReader.holeNumberArray.arrayLength, maxAllocNElements, "HoleNumbers (base)");
		holeNumbers.resize(nReads);
		zmwReader.holeNumberArray.Read(0,nReads, (unsigned int*)&holeNumbers[0]);
		return holeNumbers.size();
	}

	int Advance(int nSeq) {
		int i;
		// cannot advance past the end of this file
		if (curRead + nSeq >= nReads) { return 0; }
		for (i = curRead; i < curRead + nSeq && i < nReads; i++ ) {
 			int seqLength;
			zmwReader.numEventArray.Read(i, i+1, &seqLength);
			curBasePos += seqLength;
		}
		curRead += nSeq;
		zmwReader.Advance(nSeq);
		return curRead;
	}

	int GetNextWithoutPosAdvance(FASTASequence &seq) {
		int seqLength;

		zmwReader.numEventArray.Read(curRead, curRead+1, &seqLength);
		seq.length = 0;
		seq.seq = NULL;

		if (includedFields["Basecall"]) {
			if (seqLength > 0) {
				ResizeSequence(seq, seqLength);
				baseArray.Read(curBasePos, curBasePos + seqLength, (unsigned char*) seq.seq);
			}
		}

		string readTitle;
		unsigned int holeNumber;
		zmwReader.holeNumberArray.Read(curRead, curRead+1, &holeNumber);
		seq.StoreHoleNumber(holeNumber);

		DNALength simIndex=0, simCoordinate=0;

		if (includedFields["SimulatedSequenceIndex"] == true) {
			simulatedSequenceIndexArray.Read(curRead,curRead+1,&simIndex);
		}
		if (includedFields["SimulatedCoordinate"] == true) {
			simulatedCoordinateArray.Read(curRead, curRead+1, &simCoordinate);
		}
		
		if (scanDataReader.platformId == AstroPlatform) {
			int16_t xy[2];
			if (zmwReader.readHoleXY) {
				zmwReader.xyArray.Read(curRead, curRead+1, 0, 2, xy);
			}
			else {
				xy[0] = xy[1] = 0;
			}
			
			AstroBuildReadTitle(scanDataReader.GetRunCode(), scanDataReader.GetMovieName(), 
													xy[0], xy[1], readTitle, holeNumber, simIndex, simCoordinate);
			seq.StoreXY(xy);
		}
		else {
			SpringfieldBuildReadTitle(scanDataReader.GetMovieName(), holeNumber, readTitle, simIndex, simCoordinate);
		}
		seq.CopyTitle(readTitle);
		curRead++;
		return seqLength;
	}

	int GetNextDeletionQV(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateDeletionQVSpace(seq.length);
		deletionQVArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.deletionQV);
	}

	int GetNextDeletionTag(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateDeletionTagSpace(seq.length);
		deletionTagArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.deletionTag);
	}

	int GetNextInsertionQV(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateInsertionQVSpace(seq.length);
		insertionQVArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.insertionQV);
	}

	int GetNextWidthInFrames(SMRTSequence &seq) {
		if (seq.length == 0) return 0;
		seq.widthInFrames = new HalfWord[seq.length];
		basWidthInFramesArray.Read((int)curBasePos, (int) curBasePos + seq.length, (HalfWord*) seq.widthInFrames);
	}

	int GetNextPreBaseFrames(SMRTSequence &seq) {
		if (seq.length == 0) return 0;
		seq.preBaseFrames = new HalfWord[seq.length];
		preBaseFramesArray.Read((int)curBasePos, (int) curBasePos + seq.length, (HalfWord*) seq.preBaseFrames);
	}
	int GetNextPulseIndex(SMRTSequence &seq) {
		if (seq.length == 0) return 0;
		seq.pulseIndex = new unsigned int[seq.length];
		pulseIndexArray.Read((int)curBasePos, (int) curBasePos + seq.length, (int*) seq.pulseIndex);
	}

	int GetNextSubstitutionQV(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateSubstitutionQVSpace(seq.length);
		substitutionQVArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.substitutionQV);
	}

	int GetNextSubstitutionTag(FASTQSequence &seq) {
		if (seq.length == 0) return 0;
		seq.AllocateSubstitutionTagSpace(seq.length);
		substitutionTagArray.Read((int)curBasePos, (int) curBasePos + seq.length, (unsigned char*) seq.substitutionTag);		
	}

	void Close() {

		baseCallsGroup.Close();
		nElemArray.Close();
		zmwXCoordArray.Close();
		zmwYCoordArray.Close();
		baseArray.Close();
		qualArray.Close();
		if (useZmwReader) {
			zmwReader.Close();
		}

		if (includedFields["DeletionQV"]) {
			deletionQVArray.Close();
		}
		if (includedFields["DeletionTag"]) {
			deletionTagArray.Close();
		}
		if (includedFields["InsertionQV"]) {
			insertionQVArray.Close();
		}
		if (includedFields["SubstitutionTag"]) {
			substitutionTagArray.Close();
		}
		if (includedFields["SubstitutionQV"]) {
			substitutionQVArray.Close();
		}
		if (includedFields["WidthInFrames"]) {
			basWidthInFramesArray.Close();
		}
		if (includedFields["PreBaseFrames"]) {
			preBaseFramesArray.Close();
		}
		if (includedFields["PulseIndex"]) {
			pulseIndexArray.Close();
		}

		HDFPulseDataFile::Close();
	}

	void ReadAllHoleXY(BaseFile &baseFile) {
		baseFile.holeXY.resize(nReads);
		int i;
		for (i = 0; i < nReads; i++) {
			zmwReader.xyArray.Read(i,i+1, baseFile.holeXY[i].xy);
		}
	}

	void ReadBaseFile(BaseFile &baseFile) {
		
		if (scanDataReader.fileHasScanData) {
			scanDataReader.Read(baseFile.scanData);
		}

		baseFile.nReads = nReads;
		if (includedFields["Basecall"]) {
			baseFile.baseCalls.resize(nBases);
			baseArray.Read(0,nBases, &baseFile.baseCalls[0]);
		}

		if (useBasHoleXY) {
			ReadAllHoleXY(baseFile);
		}
		GetAllHoleNumbers(baseFile.holeNumbers);
		zmwReader.numEventArray.ReadDataset(baseFile.readLengths);
		
		/*
		 * This can probably be fixed eventually with an object factory or
		 * collection of some sorts.
		 */
		if (includedFields["WidthInFrames"]) {
			basWidthInFramesArray.ReadDataset(baseFile.basWidthInFrames);
		}
		if (includedFields["PreBaseFrames"]) {
			preBaseFramesArray.ReadDataset(baseFile.preBaseFrames);
		}
		if (includedFields["PulseIndex"]) {
			pulseIndexArray.ReadDataset(baseFile.pulseIndex);
		}
		if (includedFields["QualityValue"]) {
			qualArray.ReadDataset(baseFile.qualityValues);
		}
		if (includedFields["InsertionQV"]) {
			insertionQVArray.ReadDataset(baseFile.insertionQV);
		}
		if (includedFields["SubstitutionTag"]) {
			substitutionTagArray.ReadDataset(baseFile.substitutionTag);
		}
		if (includedFields["SubstitutionQV"]) {
			substitutionQVArray.ReadDataset(baseFile.substitutionQV);			
		}
		if (includedFields["DeletionQV"]) {
			deletionQVArray.ReadDataset(baseFile.deletionQV);
		}
		if (includedFields["DeletionTag"]) {
			deletionTagArray.ReadDataset(baseFile.deletionTag);
		}

		baseFile.nBases = nBases;
		baseFile.scanData.platformId = scanDataReader.platformId;
	}
};


template<>
int T_HDFBasReader<SMRTSequence>::Advance(int nSteps) {
	int retVal;
	retVal = ((T_HDFBasReader<FASTQSequence>*)this)->Advance(nSteps);
	return retVal;
}

typedef T_HDFBasReader<FASTASequence> HDFBasReader;
typedef T_HDFBasReader<FASTQSequence> HDFQualReader;
typedef T_HDFBasReader<SMRTSequence>  HDFSmrtReader;


#endif
