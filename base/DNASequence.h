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

#ifndef  DNA_SEQUENCE_H_
#define  DNA_SEQUENCE_H_
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <boost/shared_ptr.hpp>

#include "base/NucConversion.h"
using namespace std;
#include <iostream>
#include "base/Types.h"

class DNASequence;

typedef boost::shared_ptr<DNASequence> DNASequencePtr;

typedef u_int32_t DNALength;
typedef unsigned char Nucleotide;

class DNASequence {
 public:
  DNALength length;
  Nucleotide *seq;
	bool deleteOnExit;
	DNALength size() {
		return length;
	}
	typedef Nucleotide T_Block;
	DNASequence &Copy(const DNASequence &rhs) {
		if (length != 0) {
			assert(seq != NULL);
			delete[] seq;
		}
		seq = new Nucleotide [rhs.length];
		memcpy(seq, rhs.seq, rhs.length);
		length = rhs.length;
		deleteOnExit = true;
		return *this;
	}

	void ShallowCopy(const DNASequence &rhs) {
		seq = rhs.seq;
		length = rhs.length;
		deleteOnExit = false;
	}
	int GetStorageSize() {
		return (length * sizeof(Nucleotide));
	}

	DNASequence &operator=(const DNASequence &rhs){ 
		Copy(rhs);
		return *this;
	}

	int bitsPerNuc;
	DNASequence() {
		seq = NULL;
		length = 0;
		bitsPerNuc = 8;
		deleteOnExit = false;
	}
	//
	// synonym for printseq
	//
	void Print(ostream &out, int lineLength = 50) {
		PrintSeq(out, lineLength);
	}
	void PrintSeq(ostream &out, int lineLength = 50) {
		DNALength curPos = 0;
		int curLineLength = lineLength;
		while (curPos < length) {
			if (curPos + curLineLength > length) {
				curLineLength = length - curPos;
			}
			string line;
			line.assign((char*) &seq[curPos], curLineLength);
			out << line << std::endl;
			curPos += curLineLength;
		}
	}

	void Allocate(DNALength plength) {
		if (seq != NULL) {
			delete[] seq;
		}
		seq = new Nucleotide [plength];
		length = plength;
		deleteOnExit = true;
	}

	void ReferenceSubstring(DNASequence &rhs, int pos=0, int substrLength=0) {
		//
		// This makes a reference therefore it should not be deleted.
		//
		if (substrLength == 0) {
			substrLength = rhs.length - pos;
		}
		assert(pos + substrLength <= rhs.length);
		seq = &rhs.seq[pos];
		length = substrLength;
		deleteOnExit = false;
	}

	void MakeRC(DNASequence &rc) {
		rc.Allocate(length);
		DNALength i;
		for (i = 0; i < length; i++) {
			rc.seq[length - i - 1] = ReverseComplementNuc[seq[i]];
		}
	}

	void ToTwoBit() {
		DNALength i;
		for (i = 0; i < length; i++) {
			seq[i] = TwoBit[seq[i]];
		}
		bitsPerNuc = 2;
	}

	void ToThreeBit() {
		DNALength i;
		if (bitsPerNuc != 3) 
			for (i = 0; i < length; i++) { seq[i] = ThreeBit[seq[i]]; }
		bitsPerNuc = 3;
	}

	void ToFourBit() {
		DNALength i;
		if (bitsPerNuc != 4) 
			for (i = 0; i < length; i++) { seq[i] = FourBit[seq[i]]; }
		bitsPerNuc = 4;
	}
	
	void ConvertThreeBitToAscii() {
		DNALength i;
		for (i = 0; i < length; i++ ){
			seq[i] = ThreeBitToAscii[seq[i]];
		}
	}

	void ToAscii() {
		DNALength i;
		if (bitsPerNuc != 8) {
			for (i = 0; i < length; i++ ){ 
				seq[i] = FourBitToAscii[seq[i]];
			}
			bitsPerNuc = 8;
		}
	}
		
	void Assign(DNASequence &ref, DNALength start=0, DNALength plength=0) {
		if (seq != NULL) {
			delete[] seq;
		}
		if (plength) {
			length = plength;
			seq = new Nucleotide[length];
			memcpy(seq, &ref.seq[start], length);
		}
		else if (start) {
			length = ref.length - start;
			seq = new Nucleotide[length];
			memcpy(seq, &ref.seq[start], length);
		}
		else {
			this->Copy(ref);
		}
		deleteOnExit = true;
	}

	void ToLower() {
		DNALength i;
		for (i = 0; i < length; i++) {
			seq[i] = AllToLower[seq[i]];
		}
	}
		
	void ToUpper() {
		DNALength i;
		for (i = 0; i < length; i++) {
			seq[i] = AllToUpper[seq[i]];
		}
	}

	void Concatenate(const Nucleotide *moreSeq, DNALength moreSeqLength) {
		DNALength prevLength = length;
		length += moreSeqLength;
		Nucleotide *prev = seq;
		seq = new Nucleotide[length];
		//		cout << "adding seq " << length << endl;
		if (prev != NULL) {
			memcpy(seq, prev, prevLength);
			//			cout << "deleting prev." << endl;
			delete[] prev;
		}
		memcpy((Nucleotide*) &seq[prevLength], moreSeq, moreSeqLength);

	}

	string GetTitle() const {
		return string("");
	}

	void Concatenate(const Nucleotide* moreSeq) {
		DNALength moreSeqLength = strlen((char*) moreSeq);
		Concatenate(moreSeq, moreSeqLength);
	}
	
	void Concatenate(DNASequence &seq) {
		Concatenate(seq.seq, seq.length);
	}

	int Compare(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length) {
		return memcmp(&seq[pos], &rhs.seq[rhsPos], length);
	}

	int LessThanEqual(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length) {
		int res = Compare(pos, rhs, rhsPos, length);
		if (res <= 0) 
			return 1;
		else
			return 0;
	}

	int Equals(DNASequence &rhs, DNALength rhsPos, DNALength length, DNALength pos=0 ) {
		int res = Compare(pos, rhs, rhsPos, length);
		return res == 0;
	}

	int LessThan(DNALength pos,  DNASequence &rhs, DNALength rhsPos, DNALength length) {
		int res=  Compare(pos, rhs, rhsPos, length);
		return (res < 0);
	}

	void CleanupASCII() {
		DNALength i;
		for (i = 0; i < length; i++ ){
			if (ThreeBit[seq[i]] == 255) {
				seq[i] = 'N';
			}
		}
	}
	Nucleotide operator[](int i) {
		return GetNuc(i);
	}

	Nucleotide GetNuc(DNALength i) {
		return seq[i];
	}

	void Free() {
		if (deleteOnExit == false) { return; }
		if (seq != NULL) {
			delete[] seq;
			seq = NULL;
			length = 0;
		}
	}
	void Resize(DNALength newLength) {
		if (seq != NULL) {
			delete[] seq;
		}
		seq = new  Nucleotide[newLength];
		length = newLength;
		deleteOnExit = true;
	}
	DNALength GetSeqStorage() {
		return length;
	}
};

template<typename T>
DNALength ResizeSequence(T &dnaseq, DNALength newLength) {
	assert(newLength > 0);
	if (dnaseq.seq != NULL) {
		delete[] dnaseq.seq;
	}
	dnaseq.seq = new Nucleotide[newLength];
	dnaseq.length = newLength;
	dnaseq.deleteOnExit = true;
	return newLength;
}

#endif
