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

#ifndef NUC_CONVERSION_H_
#define NUC_CONVERSION_H_

//
// Map from ascii to 2 bit representation.
//
static int TwoBit[] = {
  	0,  1,  2,  3,  0,  1,  2,  3,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,  0,255,  1,255,255,255,  2,
  255,255,255,255,255,255,255,255,
  255,255,255,255,  3,255,255,255,
  255,255,255,255,255,255,255,255,
  255,  0,255,  1,255,255,255,  2,
  255,255,255,255,255,255,255,255,
  255,255,255,255,  3,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255
};


//
// Map from charcter to 3 bit.  Treat all non ACGT IUPAC as N
// IPUAC: ACGT U M R W S Y K V H D N
// 

static int ThreeBit[] = {
  	0,  1,  2,  3,    4,255,255,255, // 0..3 are native nucletides, 4
																		 // is an 'N', which is
																		 // represented just fine in the 3
																		 // bit mode.
		255,255,255,255,255,255,255,255, // 8
		255,255,255,255,255,255,255,255, // 16
		255,255,255,255,255,255,255,255, // 24
		255,255,255,255,  5,255,255,255, // 32  define '$' for bwt.
		255,255,255,255,255,255,255,255, // 40
		255,255,255,255,255,255,255,255, // 48
		255,255,255,255,255,255,255,255, // 56
		255,  0,  4,  1,  4,255,255,  2, // 64
		  4,255,255,  4,255,  4,  4,255, // 72
		255,255,  4,  4,  3,  4,  4,  4, // 80
		255,  4,255,255,255,255,255,  4, // 88
		255,  0,  4,  1,  4,255,255,  2, // 96
		  4,255,255,  4,255,  4,  4,255, // 104
		255,255,  4,  4,  3,  4,  4,  4, // 112
		  4,  4,255,255,255,255,255,255, // 120
		255,255,255,255,255,255,255,255, // 128
		255,255,255,255,255,255,255,255, // 136
		255,255,255,255,255,255,255,255, // 144
		255,255,255,255,255,255,255,255, // 152
		255,255,255,255,255,255,255,255, // 160
		255,255,255,255,255,255,255,255, // 168
		255,255,255,255,255,255,255,255, // 176
		255,255,255,255,255,255,255,255, // 184
		255,255,255,255,255,255,255,255, // 192
		255,255,255,255,255,255,255,255, // 200
		255,255,255,255,255,255,255,255, // 208
		255,255,255,255,255,255,255,255, // 216
		255,255,255,255,255,255,255,255, // 224
		255,255,255,255,255,255,255,255, // 232
		255,255,255,255,255,255,255,255, // 240
		255,255,255,255,255,255,255,255  // 248
};

static int IsACTG[] = {
   	1,  1,  1,  1,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  1,  0,  1,  0,  0,  0,  1,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  1,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  1,  0,  1,  0,  0,  0,  1,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  1,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0
};


static int FourBit[] = {
   	0,  1,  2,  3,  4,  5,  6,  7,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,  0,255,  1,255,255,255,  2,
  255,255,255,255,255,255,  8,255,
  255,255,255,255,  3,255,255,255,
  255,255,255,255,255,255,255,255,
  255,  4,255,  5,255,255,255,  6,
  255,255,255,255,255,255,  8,255,
  255,255,255,255,  7,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255
};

static int FourBitToAscii[]  = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'N', 'X'};
static int ThreeBitToAscii[] = {'A', 'C', 'G', 'T', 'N', '$'};

static int MaskedFourBit[] = {
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,0,9,1,9,9,9,2,
  9,9,9,9,9,9,8,9,
  9,9,9,9,3,9,9,9,
  9,9,9,9,9,9,9,9,
  9,4,9,5,9,9,9,6,
  9,9,9,9,9,9,8,9,
  9,9,9,9,7,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9,
  9,9,9,9,9,9,9,9
};

static int AllToUpper[] = {
	'A','C','G','T','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','A','N','C','N','N','N','G',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','T','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','A','N','C','N','N','N','G',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','T','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N'
};
static int AllToLower[] = {
	'a','c','g','t','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','a','N','c','N','N','N','g',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','t','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','a','N','c','N','N','N','g',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','t','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N'
};


static char TwoBitToAscii[] = {'A', 'C', 'G', 'T'};


static unsigned char ReverseComplementNuc[] = {
    3,  2,  1,  0,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,'T',127,'G',127,127,127,'C',
	127,127,127,127,127,127,'N',127,
  127,127,127,127,'A',127,127,127,
  127,127,127,127,127,127,127,127,
  127,'t',127,'g',127,127,127,'c',
  127,127,127,127,127,127,'n',127,
  127,127,127,127,'a',127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127
};

static int NucToHdfColumnOrder[] = {
  	2,  3,  1,  0,  2,  3,  1,  0,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,  2,255,  3,255,255,255,  1,
  255,255,255,255,255,255,255,255,
  255,255,255,255,  0,255,255,255,
  255,255,255,255,255,255,255,255,
  255,  2,255,  3,255,255,255,  1,
  255,255,255,255,255,255,255,255,
  255,255,255,255,  0,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255,
  255,255,255,255,255,255,255,255
};

#endif
