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

#ifndef DATA_HDF_HDF_REGION_TABLE_READER_H_
#define DATA_HDF_HDF_REGION_TABLE_READER_H_

#include "datastructures/reads/RegionTable.h"
#include "data/hdf/HDFFile.h"
#include "data/hdf/HDFArray.h"
#include "data/hdf/HDF2DArray.h"
#include "data/hdf/HDFAtom.h"
#include "data/hdf/PlatformId.h"
#include <string>

using namespace H5;
using namespace std;


class HDFRegionTableReader {
 public:
	H5File regionTableFile;
	HDFGroup pulseDataGroup;
	HDF2DArray<int> regions;
	
	HDFAtom<vector<string> > regionTypes;
	HDFAtom<vector<string> > regionDescriptions;
	HDFAtom<vector<string> > regionSources;
	HDFAtom<vector<string> > columnNames;
	int curRow;
	int nRows;
	bool fileContainsRegionTable;

	int Initialize(string &regionTableFileName) {
		/*
		 * Initialize access to the HDF file.
		 */
		try {
			regionTableFile.openFile(regionTableFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		}
		catch (Exception &e) {
			cout << e.getDetailMsg() << endl;
			return 0;
		}
		if (pulseDataGroup.Initialize(regionTableFile, "PulseData") == 0) {
			return 0;
		}
		if (pulseDataGroup.ContainsObject("Regions") == 0) {
			fileContainsRegionTable = false;
			return 0;
		}
		else {
			fileContainsRegionTable = true;
		}

		if (regions.Initialize(regionTableFile, "PulseData/Regions") == 0) {
			return 0;
		}
		nRows = regions.GetNRows();

		if (columnNames.Initialize(regions.dataset, "ColumnNames") == 0) {
			return 0;
		}
		if (regionTypes.Initialize(regions.dataset, "RegionTypes") == 0) {
			return 0;
		}
		if (regionDescriptions.Initialize(regions.dataset, "RegionDescriptions") == 0) {
			return 0;
		}
		if (regionSources.Initialize(regions.dataset,  "RegionSources") == 0) {
			return 0;
		}
		
		curRow = 0;
		return 1;
	}

	int GetNext(RegionAnnotation &annotation) {
		//
		// Bail with no-op if this is the last row.
		//

		if (fileContainsRegionTable == false) {
			return 0;
		}

		if (curRow == nRows) {
			return 0;
		}

		regions.Read(curRow, curRow+1, annotation.row);
		++curRow;
		return 1;
	}	

	void RegionTypesToMap(RegionTable &table) {
		int i;
		table.regionTypeEnums.resize(table.regionTypes.size());
		for (i = 0;i < table.regionTypes.size(); i++) {
			if (table.regionTypes[i] == "GlobalAccuracy") {
				table.regionTypeEnums[i] = GlobalAccuracy;
			}
			else if (table.regionTypes[i] == "HQRegion") {
				table.regionTypeEnums[i] = HQRegion;
			}
			else if (table.regionTypes[i] == "Adapter") {
				table.regionTypeEnums[i] = Adapter;
			}
			else if (table.regionTypes[i] == "Insert") {
				table.regionTypeEnums[i] = Insert;
			}
			else if (table.regionTypes[i] == "Accuracy") {
				table.regionTypeEnums[i] = Insert;
			}
			else {
				cout << "ERROR! Region Type " << table.regionTypes[i] << " is not supported.  Check Enumerations.h" << endl;
				assert(0);
			}
		}
	}

	int ReadTableAttributes(RegionTable &table) {
		if (fileContainsRegionTable == false) {
			return 0;
		}
		columnNames.Read(table.columnNames);
		regionTypes.Read(table.regionTypes);
		RegionTypesToMap(table);
		regionDescriptions.Read(table.regionDescriptions);
		regionSources.Read(table.regionSources);
	}
	
	void Close() {
		pulseDataGroup.Close();
		regions.Close();
		regionTableFile.close();
	}

	void ReadTable(RegionTable &table) {
		if (fileContainsRegionTable == false) {
			return;
		}
		ReadTableAttributes(table);
		table.table.resize(nRows);
		int i;
		while(GetNext(table.table[curRow])) {
			i++;
		}
	}
	

};


#endif
