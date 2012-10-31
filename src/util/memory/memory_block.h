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
#ifndef MEMORY_BLOCK_H
#define MEMORY_BLOCK_H

#include <iterator>
#include "common.h"


namespace util {

struct MemoryBlock {
    uint length;
    byte * mem;
    MemoryBlock( byte * m, uint len) : mem(m), length( len ) {}
};

/**
    MemoryBlockIterator is an abstract iterator for accessing sub-blocks within a MemoryBlock

    The idea being that a memory block can contain
**/
class MemoryBlockIterator : public std::iterator< input_iterator_tag, byte > {
    public:
        MemoryBlockIterator( MemoryBlock * mb ) : ptr( mb->mem ), length( mb->length ), cur_idx(0), sub_block_size(0) {}
        MemoryBlockIterator( const MemoryBlockIterator& mbi ) : ptr( mbi.ptr ), length( mbi.length ), cur_idx( mbi.cur_idx ), sub_block_size( mbi.sub_block_size ) {}

        MemoryBlockIterator& operator++() { next(); return *this; }

        inline bool operator==( const MemoryBlockIterator& mbi ) { return ptr == mbi.ptr && sub_block_size == mbi.sub_block_size; }
        bool operator!=( const MemoryBlockIterator& mbi ) { return !operator==(mbi); }

        MemoryBlock& operator*() { return MemoryBlock( ptr, sub_block_size ); }

    protected:
        virtual void next() = 0;
        int length, cur_idx, sub_block_size;
        byte * ptr;
};

}

#endif // MEMORY_BLOCK_H
