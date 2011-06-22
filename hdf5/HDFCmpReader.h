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

#ifndef DATA_HDF_HDF_CMP_READER_H_
#define DATA_HDF_HDF_CMP_READER_H_

#include "H5Cpp.h"
#include <iostream>
#include <assert.h>
#include "datastructures/alignment/CmpFile.h"
#include "datastructures/alignment/Alignment.h"
#include "datastructures/alignment/CmpAlignment.h"
#include "datastructures/alignment/CmpReadGroupTable.h"
#include "datastructures/alignment/CmpRefSeqTable.h"
#include "data/hdf/HDFAtom.h"
#include "data/hdf/HDFArray.h"
#include "data/hdf/HDF2DArray.h"
#include "data/hdf/HDFCmpData.h"
#include "data/hdf/HDFCmpRefAlignmentGroup.h"
#include "data/hdf/HDFCmpExperimentGroup.h"
#include "HDFAlnGroupGroup.h"
#include "HDFAlnInfoGroup.h"
#include "HDFRefGroupGroup.h"
#include "HDFRefInfoGroup.h"
#include "HDFMovieInfoGroup.h"
#include "HDFCmpRootGroup.h"
#include <sstream>
#include <map>

using namespace H5;
using namespace std;

template <typename T_Alignment>
class HDFCmpReader : public HDFCmpData {
public:
	map<int,int>  movieNameIdToArrayIndex,  readGroupPathIdToArrayIndex, refGroupIdToArrayIndex;
	map<string,string> readGroupPathToReadGroup;
	HDF2DArray<int> alignmentIndexArray;
	HDFAlnGroupGroup alnGroupGroup;
	HDFAlnInfoGroup  alnInfoGroup;
	HDFMovieInfoGroup movieInfoGroup;
	HDFRefGroupGroup refGroupGroup;
	HDFRefInfoGroup refInfoGroup;
	HDFCmpRootGroup<T_Alignment> rootGroup;

	void AstroInitializeColumnNameMap() {
		CmpAlignment::columnNameToIndex["AlignmentId"] = 0;
		CmpAlignment::columnNameToIndex["ReadGroupId"] = 1;
		CmpAlignment::columnNameToIndex["MovieId"] = 2;
		CmpAlignment::columnNameToIndex["RefSeqId"] = 3;		
		CmpAlignment::columnNameToIndex["tStart"] = 4;
		CmpAlignment::columnNameToIndex["tEnd"] = 5;
		CmpAlignment::columnNameToIndex["AlignedStrand"] = 6;
		CmpAlignment::columnNameToIndex["ExpId"] = 7;
		CmpAlignment::columnNameToIndex["RunId"] = 8;		
		CmpAlignment::columnNameToIndex["Panel"] = 9;
		CmpAlignment::columnNameToIndex["x"] = 10;
		CmpAlignment::columnNameToIndex["y"] = 11;
		CmpAlignment::columnNameToIndex["SubreadId"] = 12;
		CmpAlignment::columnNameToIndex["rStart"] = 13;		
		CmpAlignment::columnNameToIndex["rEnd"] = 14;
		CmpAlignment::columnNameToIndex["Z"] = 15;
		CmpAlignment::columnNameToIndex["nM"] = 16;
		CmpAlignment::columnNameToIndex["nMM"] = 17;
		CmpAlignment::columnNameToIndex["nIns"] = 18;		
		CmpAlignment::columnNameToIndex["nDel"] = 19;
		CmpAlignment::columnNameToIndex["offset_begin"] = 20;
		CmpAlignment::columnNameToIndex["offset_end"] = 21;
	}

	void SpringfieldInitializeColumnNameMap() {
		CmpAlignment::columnNameToIndex["AlignmentId"] = 0;
		CmpAlignment::columnNameToIndex["ReadGroupId"] = 1;
		CmpAlignment::columnNameToIndex["MovieId"] = 2;
		CmpAlignment::columnNameToIndex["RefSeqId"] = 3;		
		CmpAlignment::columnNameToIndex["tStart"] = 4;
		CmpAlignment::columnNameToIndex["tEnd"] = 5;
		CmpAlignment::columnNameToIndex["RCRefStrand"] = 6;
		CmpAlignment::columnNameToIndex["HoleNumber"] = 7;
		CmpAlignment::columnNameToIndex["SetNumber"] = 8;		
		CmpAlignment::columnNameToIndex["StrobeNumber"] = 9;
		CmpAlignment::columnNameToIndex["SubreadId"] = 10;
		CmpAlignment::columnNameToIndex["rStart"] = 11;		
		CmpAlignment::columnNameToIndex["rEnd"] = 12;
		CmpAlignment::columnNameToIndex["MapQV"] = 13;
		CmpAlignment::columnNameToIndex["nM"] = 14;
		CmpAlignment::columnNameToIndex["nMM"] = 15;
		CmpAlignment::columnNameToIndex["nIns"] = 16;		
		CmpAlignment::columnNameToIndex["nDel"] = 17;
		CmpAlignment::columnNameToIndex["offset_begin"] = 18;
		CmpAlignment::columnNameToIndex["offset_end"] = 19;
		CmpAlignment::columnNameToIndex["nBackRead"] = 20;
		CmpAlignment::columnNameToIndex["nReadOverlap"] = 21;
	}


	int Initialize(string &hdfCmpFileName, unsigned int flags=H5F_ACC_RDONLY) {
		/*
		 * Initialize access to the HDF file.
		 */
		try {
			hdfCmpFile.openFile(hdfCmpFileName.c_str(), flags, H5P_DEFAULT);
		}
		catch (Exception &e) {
			cout << e.getDetailMsg() << endl;
			return 0;
		}

		rootGroup.Initialize(hdfCmpFile);
		
		if (alnGroupGroup.Initialize(rootGroup.rootGroup) == 0)       { return 0; }
		if (refInfoGroup.Initialize(rootGroup.rootGroup) == 0)   { return 0; }
		if (refGroupGroup.Initialize(rootGroup.rootGroup) == 0)   { return 0; }
		if (movieInfoGroup.Initialize(rootGroup.rootGroup) == 0) { return 0; }
		if (alnInfoGroup.Initialize(rootGroup.rootGroup) == 0)   { return 0; }

		return 1;
	}

	static void ParseReadGroupPath(string &path, string &refName, string &readGroupName){ 
		int delimPos;
		delimPos = path.find_last_of('/');
		if (delimPos != path.npos) {
			refName = path.substr(0, delimPos);
			//
			// Check the ref name to see if it has a preceeding '/'
			//
			int firstSlashPos;
			firstSlashPos = refName.find_first_of('/');
			readGroupName = path.substr(delimPos+1, path.size());
		} 
		else {
			refName ="";
			readGroupName = "";
		}
	}

	void StorePlatformId(CmpFile<T_Alignment> &cmpFile) {
		if (cmpFile.colNames[6] == "AlignedStrand") {
			cmpFile.platformId = Astro;
		}
		else {
			cmpFile.platformId = Springfield;
		}
	}


	void Read(CmpFile<T_Alignment> &cmpFile) {
		
		//
		// Gather run information.
		//
		rootGroup.ReadAttributes(cmpFile);

		//
		// Read in the groups that describe what alignments are present.
		//
		
		alnGroupGroup.Read(cmpFile.alnGroup);
		alnInfoGroup.Read(cmpFile.alnInfo);
		refGroupGroup.Read(cmpFile.refGroup);
		movieInfoGroup.Read(cmpFile.movieInfo);
		refInfoGroup.Read(cmpFile.refInfo);

		/*
		 * Now for every reference group in the cmp file, create a group.
		 */

		map<string,int> refNameToArrayIndex;

		unsigned int refSeqIndex;

		for (refSeqIndex = 0; refSeqIndex < cmpFile.refGroup.path.size(); refSeqIndex++) {
			HDFCmpRefAlignmentGroup* refAlignGroup;
			refAlignGroup = new HDFCmpRefAlignmentGroup;
			refAlignGroup->Initialize(rootGroup.rootGroup.group, cmpFile.refGroup.path[refSeqIndex]);
			int refAlignGroupIndex = refAlignGroups.size();
			refAlignGroups.push_back(refAlignGroup);				
			//
			// Allow one to index directly into this group given a string
			// representing the reference name.
			//
			refNameToArrayIndex[cmpFile.refGroup.path[refSeqIndex]] = refAlignGroupIndex;
			refGroupIdToArrayIndex[cmpFile.refGroup.id[refSeqIndex]] = refAlignGroupIndex;
		}

		//
		// Build a map from movie name to index.  This allows translation
		// of the group name to an index when storing all movie groups
		// that exist under ref groups.
		//
		int movieIndex;
		map<string,int> movieNameToId;
		for (movieIndex = 0; movieIndex < cmpFile.movieInfo.name.size(); movieIndex++) {
			movieNameToId[cmpFile.movieInfo.name[movieIndex]] = cmpFile.movieInfo.id[movieIndex];
		}

		/*
		 * Now that all ref groups are created, add groups for movies that
		 * are aligned to the ref groups.
		 */

		unsigned int readGroupIndex;
		vector<unsigned char> alignmentArray;
		for (readGroupIndex = 0; readGroupIndex < cmpFile.alnGroup.path.size(); readGroupIndex++) {
			string refName, movieName;
			ParseReadGroupPath(cmpFile.alnGroup.path[readGroupIndex], refName, movieName);

			//
			// Create an index that allows one to immediately find the movie
			// name given a readGroupIndex
			//
			readGroupPathToReadGroup[cmpFile.alnGroup.path[readGroupIndex]] = movieName;

			//
			// Look up the group where this alignment should be found.
			unsigned int refGroupArrayIndex;
			if (refNameToArrayIndex.find(refName) != refNameToArrayIndex.end()) {
				refGroupArrayIndex = refNameToArrayIndex[refName];
			}
			else {
				cout << "ERROR! The reference name '" << refName << "' does not have an entry though it is "
						 << " specified in the path for " << cmpFile.readGroupTable.names[readGroupIndex] 
						 << endl;
				assert(0);
			}
			HDFCmpExperimentGroup *alnGroupPtr;
			alnGroupPtr = refAlignGroups[refGroupArrayIndex]->InitializeExperimentGroup(movieName);
			refAlignGroups[refGroupArrayIndex]->movieIdToIndex[movieNameToId[movieName]] = refAlignGroups[refGroupArrayIndex]->readGroups.size()-1;
		}
		
		/*
		 * Now that the alignment indices are all read in, read the base-by-base alignments.
		 */
		
		
		unsigned int alignmentIndex;

		for (alignmentIndex = 0; alignmentIndex < cmpFile.alnInfo.alignments.size(); alignmentIndex++) {
			int alnGroupId = cmpFile.alnInfo.alignments[alignmentIndex].GetAlnGroupId();
			int refGroupId = cmpFile.alnInfo.alignments[alignmentIndex].GetRefGroupId();
			int movieId    = cmpFile.alnInfo.alignments[alignmentIndex].GetMovieId();
			string refSeqName, readGroupName;

			//
			// Make sure the refGroupId specified in the alignment index exists in the alignment groups.
			//
			int refGroupArrayIndex;
			if (refGroupIdToArrayIndex.find(refGroupId) == refGroupIdToArrayIndex.end()) {
					cout << "ERROR! Alignment " << cmpFile.alnInfo.alignments[alignmentIndex].GetAlignmentId()
							 << " has ref seq id " << refGroupId << " that does not exist in the HDF file." << endl;
					assert(0);
			}
			else {
				refGroupArrayIndex = refGroupIdToArrayIndex[refGroupId];
			}
			
			//
			// Point to the refGroup that this alignment is part of.
			//
			HDFCmpRefAlignmentGroup* refAlignGroup = refAlignGroups[refGroupArrayIndex];

			//
			// Now locate the movie group that is part of this ref align group.
			//
			string movieName;
			if (cmpFile.movieInfo.FindMovie(movieId, movieName) == 0) {
				cout << "ERROR! Alignment " << cmpFile.alnInfo.alignments[alignmentIndex].GetAlignmentId() 
						 << " specifies a movie id " << movieId << " that is not listed in the movie group." << endl;
				assert(0);
			}
			
			if (refAlignGroup->experimentNameToIndex.find(movieName) ==
					refAlignGroup->experimentNameToIndex.end()) {
				cout << "Internal ERROR! The movie name " << movieName << " is specified as part of "
						 << " the path in alignment " << cmpFile.alnInfo.alignments[alignmentIndex].GetAlignmentId() 
						 << " though it does not exist in the ref align group specified for this alignment." << endl;
				assert(0);
			}
			
			int experimentIndex;
			experimentIndex = refAlignGroup->movieIdToIndex[movieId];

			int alignmentLength = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd() - 
				cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
			if (alignmentArray.size() < alignmentLength) {
				alignmentArray.resize(alignmentLength);
			}
 			refAlignGroup->readGroups[experimentIndex]->alignmentArray.Read(cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin(),
																																			cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd(),
																																			&alignmentArray[0]);

			cmpFile.alnInfo.alignments[alignmentIndex].StoreAlignmentArray(&alignmentArray[0], alignmentLength);
		}

	}


	void ReadReadGroupPathTable(CmpIndexedStringTable &readGroupPathTable) {
		readGroupPathIdLastRow.Read(readGroupPathTable.lastRow);
		readGroupPathTable.resize(readGroupPathTable.lastRow); // resizes all tables
		readGroupPathIdArray.Read(0, readGroupPathTable.lastRow, &readGroupPathTable.ids[0]);
		readGroupPathArray.Read(0, readGroupPathTable.lastRow, &readGroupPathTable.names[0]);
		readGroupPathTable.StoreArrayIndexMap();
	}

	void ReadRefSeqTable(CmpIndexedStringTable &refSeqTable) {
		refSeqNameIdLastRow.Read(refSeqTable.lastRow);
		refSeqTable.resize(refSeqTable.lastRow);
		refSeqNameIdArray.Read(0, refSeqTable.lastRow, &refSeqTable.ids[0]);
		refSeqNameArray.Read(0, refSeqTable.lastRow, &refSeqTable.names[0]);
		refSeqTable.StoreArrayIndexMap();
	}

	void ReadMovieNameTable(CmpIndexedStringTable &movieTable) {
		movieNameIdLastRow.Read(movieTable.lastRow);
		movieTable.resize(movieTable.lastRow);
		movieNameIdArray.Read(0, movieTable.lastRow, &movieTable.ids[0]);
		movieNameArray.Read(0,   movieTable.lastRow, &movieTable.names[0]);
		movieTable.StoreArrayIndexMap();
	}

	
	void ReadAlignmentIndices(CmpFile<T_Alignment> &cmpFile) {
	}
};


#endif
