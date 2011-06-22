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

#ifndef DATA_HDF_HDF_RAW_PULSE_READER_H_
#define DATA_HDF_HDF_RAW_PULSE_READER_H_

//
// Load pulse information from the 
//DEFAULT_PULSE_FIELDS = [ "PulseWidth", "pkmid", "StartTime", "IPD", "Light", "ClassifierQV", "QualityValue" ]
//                                     ,        ,  u32       ,      
//BASEMAP = numpy.array(['-','A','C','-','G','-','-','-','T', '-','-','-','-','-','-','-','N'])


#include "data/hdf/HDFArray.h"
#include "data/hdf/HDFFile.h"

class HDFRawPulseReader {
 public:
	vector<HDFArray*> fields;
	vector<string> fieldNames;
	H5File hdfRawPulseFile;
	Group  rootGroup;
	Initialize(string fileName) {
		/*
		 * Open the file for reading.
		 */
		try {
			hdfRawPulseFile.openFile(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		}
		catch (Exception &e) {
			cout << e.getDetailMsg() << endl;
			return 0;
		}
		rootGroup = hdfRawPulseFile.openGroup("/");
		
		fieldNames.push_back("PulseWidth");
		fieldNames.push_back("pkmid");
		fieldNames.push_back("StartTime");
		fieldNames.push_back("IPD");
		fieldNames.push_back("Light");
		fieldNames.push_back("ClassifierQV");
		fieldNames.push_back("QualityValue");
		
		int filedIndex;
		for (fieldIndex = 0; fieldIndex < fieldNames.size(); fieldIndex++) {
			fields.push_back(new HDFArray);
		}
	}
		


};


#endif


