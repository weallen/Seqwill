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
  kChrX = 21,
  kChrY = 22,
  kChrM = 23,
  kChrUnknown = 24
};

static const int NUM_CHRS = 21;

ChromosomeEnum chroms[NUM_CHRS] = {kChr1, 
  kChr2, kChr3, kChr4, kChr5, kChr6, kChr7, kChr8,                                   
  kChr9, kChr10, kChr11, kChr12, kChr13, kChr14, kChr15,
  kChr16, kChr17, kChr18, kChrX, kChrY, kChrM};

ChromosomeEnum 
ChrToNum(const std::string& chr) 
{
  if (chr.compare("chr1") == 0 || chr.compare("Chr1") == 0) 
    return kChr1;
  if (chr.compare("chr2")  == 0 || chr.compare("Chr2") == 0) 
    return kChr2;
  if (chr.compare("chr3")  == 0|| chr.compare("Chr3") == 0) 
    return kChr3;
  if (chr.compare("chr4") == 0 || chr.compare("Chr4")== 0) 
    return kChr4;
  if (chr.compare("chr5") == 0 || chr.compare("Chr5") == 0) 
    return kChr5;
  if (chr.compare("chr6") == 0 || chr.compare("Chr6") == 0) 
    return kChr6;
  if (chr.compare("chr7") == 0 || chr.compare("Chr7") == 0) 
    return kChr7;
  if (chr.compare("chr8") == 0 || chr.compare("Chr8") == 0) 
    return kChr8;
  if (chr.compare("chr9") == 0 || chr.compare("Chr9") == 0) 
    return kChr9;
  if (chr.compare("chr10") == 0 || chr.compare("Chr10") == 0) 
    return kChr10;
  if (chr.compare("chr11") == 0 || chr.compare("Chr11") == 0) 
    return kChr11;
  if (chr.compare("chr12") == 0 || chr.compare("Chr12") == 0) 
    return kChr12;
  if (chr.compare("chr13") == 0 || chr.compare("Chr13") == 0) 
    return kChr13;
  if (chr.compare("chr14") == 0 || chr.compare("Chr14") == 0) 
    return kChr14;
  if (chr.compare("chr15") == 0 || chr.compare("Chr15") == 0) 
    return kChr15;  
  if (chr.compare("chr16") == 0 || chr.compare("Chr16") == 0) 
    return kChr16;
  if (chr.compare("chr17") == 0 || chr.compare("Chr17") == 0) 
    return kChr17;
  if (chr.compare("chr18") == 0 || chr.compare("Chr18") == 0) 
    return kChr18;
  if (chr.compare("chrX") == 0 || chr.compare("ChrX") == 0) 
    return kChrX;
  if (chr.compare("chrY") == 0 || chr.compare("ChrY") == 0) 
    return kChrY;
  if (chr.compare("chrM") == 0 || chr.compare("ChrM") == 0) 
    return kChrM;
  return kChrUnknown;
} 
/*
std::string ChrToString(ChromosomeEnum chr) 
{
  if (chr == kChr1)
    return "Chr1";
  if (chr == kChr2)
    return "Chr2";
  if (chr == kChr3)
    return "Chr3";
  if (chr == kChr4)
    return "Chr4";
  if (chr == kChr5)
    return "Chr5";
  if (chr == kChr6)
    return "Chr6";
  if (chr == kChr7)
    return "Chr7";
  if (chr == kChr8)
    return "Chr8";
  if (chr == kChr9)
    return "Chr9";
  if (chr == kChr10)
    return "Chr10";
  if (chr == kChr11)
    return "Chr11";
  if (chr == kChr12)
    return "Chr12";
  if (chr == kChr13)
    return "Chr13";
  if (chr == kChr14)
    return "Chr14";
  if (chr == kChr15)
    return "Chr15";
  if (chr == kChr16)
    return "Chr16";
  if (chr == kChr17)
    return "Chr17";
  if (chr == kChr18)
    return "Chr18";
  if (chr == kChrX)
    return "ChrX";
  if (chr == kChrY)
    return "ChrY";
  if (chr == kChrM)
    return "ChrM";
  else
    return "Unknown";
}
*/
#endif
