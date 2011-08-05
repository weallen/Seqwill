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

#ifndef TYPES_H_
#define TYPES_H_

#include <sys/types.h>
#include <string>

typedef unsigned long ULong;

typedef struct PlusMinusDataFloat
{
  float plus;
  float minus;
}; 

typedef struct PlusMinusDataInt
{
  int plus;
  int minus;
};

//
// Add definitions to handle 64/32 bit computing environments
//

typedef u_int32_t VectorIndex;
typedef u_int32_t UInt;
typedef u_int8_t  Byte;
typedef u_int16_t HalfWord;

enum ChromosomeEnum {
  kChr1 = 1,
  kChr2 = 2,
  kChr3 = 3,
  kChr4 = 4,
  kChr5 = 5,
  kChr6 = 6,
  kChr7 = 7,
  kChr8 = 8,
  kChr9 = 9,
  kChr10 = 10,
  kChr11 = 11,
  kChr12 = 12,
  kChr13 = 13,
  kChr14 = 14,
  kChr15 = 15,
  kChr16 = 16,
  kChr17 = 17,
  kChr18 = 18,
  kChr19 = 19,
  kChrX = 21,
  kChrY = 22,
  kChrM = 23,
  kChrUnknown = 24
};

static const int NUM_CHRS = 22;

#endif
