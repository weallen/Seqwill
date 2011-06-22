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

#ifndef DATA_HDF_HDF_ZMW_READER_H_
#define DATA_HDF_HDF_ZMW_READER_H_

#include "datastructures/reads/ZMWGroupEntry.h"
#include "data/hdf/HDFGroup.h"

class HDFZMWReader {
 public:
	HDFGroup *parentGroupPtr;
	HDFGroup zmwGroup;
	HDFArray<unsigned int> holeNumberArray;
	HDF2DArray<int16_t> xyArray;
	HDFArray<int> numEventArray;
	int curRead;
	bool readHoleNumber;
	bool readHoleXY;
	bool readNumEvent;
	int  curZMW;
	int  nZMWEntries;
	bool  closeFileOnExit;
	H5File hdfPlsFile;
	HDFZMWReader() {
		closeFileOnExit = false;
		readHoleNumber = false;
		readHoleXY          = false;
		readNumEvent    = false;
	}

	int Initialize(HDFGroup *parentGroupP) {
		parentGroupPtr = parentGroupP;
		closeFileOnExit = false;
		return Initialize();
	}
	/*
	 * Remove this functionality for now -- assume that zmw readers must
	 * exist inside other readers such as base or pls readers.

	int InitializeFromFile(string fileName, string pathToZMWGroup) {

		destroyReaderOnExit = true;
		
		try {
			hdfPlsFile.openFile(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		}
		catch (Exception &e) {
			cout << e.getDetailMsg() << endl;
			return 0;
		}
		if (parentGroup.Initialize(hdfPlsFile, pathToZMWGroup) == 0) {
			return 0;
		}
		
		Initialize(pathToZMWGroup);
	}
	*/
	int Initialize() {
		curRead = 0;
		
		//
		// Make sure we can open the component containing the zmw information.
		//
		if (parentGroupPtr->ContainsObject("ZMW") == 0 or
				zmwGroup.Initialize(parentGroupPtr->group, "ZMW") == 0) {
			return 0;
		}
		
		//
		// Now open all the important datasets in the zmw group.  Some of
		// these are optional, so flags must be set if they do not exist.
		//
		if (zmwGroup.ContainsObject("HoleNumber")) {
			if (holeNumberArray.Initialize(zmwGroup, "HoleNumber") == 0) {
				return 0;
			}
			readHoleNumber = true;
		}
		else {
			readHoleNumber = false;
		}
		if (zmwGroup.ContainsObject("HoleXY")) {
			if (xyArray.Initialize(zmwGroup, "HoleXY") == 0) {
				return 0;
			}
			readHoleXY = true;
		}
		else {
			readHoleXY = false;
		}
		if (numEventArray.Initialize(zmwGroup, "NumEvent") == 0) {
			return 0;
		}
		nZMWEntries = numEventArray.arrayLength;
		readNumEvent = true;
		curZMW      = 0;
		return 1;
	}

	int Advance(int nSteps) {
		if (curZMW >= nZMWEntries) {
			return 0;
		}
		else {
			curZMW += nSteps;
			return 1;
		}
	}

	bool GetNext(ZMWGroupEntry &groupEntry) {
		if (curZMW == nZMWEntries) {
			return false;
		}
		if (readHoleNumber) {
			holeNumberArray.Read(curZMW, curZMW+1, &groupEntry.holeNumber);
		}
		if (readHoleXY){ 
			int16_t holeXY[2];
			xyArray.Read(curZMW, curZMW+1, holeXY);
			groupEntry.x = holeXY[0];
			groupEntry.y = holeXY[1];
		}
		numEventArray.Read(curZMW, curZMW+1, &groupEntry.numEvents);
		curZMW++;
		return true;
	}

	void Close() {
		if (readHoleNumber) {
			holeNumberArray.Close();
		}
		if (readHoleXY) {
			xyArray.Close();
		}
		if (readNumEvent) {
			numEventArray.Close();
		}

		if (closeFileOnExit == true) {
			//
			// This instance is owner of it's reader.  Close the reader file.
			//
			hdfPlsFile.close();
		}
		zmwGroup.Close();

	}
	~HDFZMWReader() {
		Close();
	}

};

#endif
