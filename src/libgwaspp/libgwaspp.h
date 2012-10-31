/*
* Copyright (c) 2012, Patrick Putnam
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met: 
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer. 
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution. 
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The views and conclusions contained in the software and documentation are those
* of the authors and should not be interpreted as representing official policies, 
* either expressed or implied, of the FreeBSD Project.
*/
#ifndef BIOLIB_H
#define BIOLIB_H

#include "libgwasppConfig.h"

#include "common.h"

#include <iostream>
#include <cstring>
#include <fstream>

#if PROCESSOR_WORD_SIZE == 64
#define PWORD ulong
#elif PROCESSOR_WORD_SIZE == 32
#define PWORD uint
#elif PROCESSOR_WORD_SIZE == 16
#define PWORD ushort
#else
#error "Expected use on 32 or 64 bit processors"
#endif

struct StringPtrComparer {
    bool operator() (const std::string * lhs, const std::string * rhs) const { return (lhs != rhs) && std::strcmp( lhs->c_str(), rhs->c_str()) < 0; }
    bool operator() (const std::string & lhs, const std::string & rhs) const { return std::strcmp( lhs.c_str(), rhs.c_str()) < 0; }
};

#endif // VERSION_H
