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
#ifndef INDEXER_H
#define INDEXER_H

#include <iostream>
#include <vector>

#include "util/index_set/index_iterator.h"

namespace util {

class indexer {
    public:

        virtual int orderOf( const std::string & id ) = 0;
        virtual int indexOf( int ord ) = 0;

        virtual std::string getIDAtOrderedIndex( int ord ) = 0;

        virtual void include( int idx, int at ) = 0;
        virtual void include( const std::vector<int>& idx ) = 0;
        virtual void include( int *indices, int count ) = 0 ;

        virtual void include( const std::string &id, int at ) = 0;
        virtual void include( const std::vector<std::string>& id ) = 0;

        virtual void exclude( int idx ) = 0;
        virtual void exclude( const std::vector<int>& idx ) = 0;
        virtual void exclude( int *indices, int count ) = 0;

        virtual void exclude( const std::string &id ) = 0;
        virtual void exclude( const std::vector<std::string>& id ) = 0;

        virtual IndexIterator *included_begin()  = 0;
        virtual IndexIterator *included_end()  = 0;

        virtual IndexIterator *excluded_begin()  = 0;
        virtual IndexIterator *excluded_end()  = 0;

        virtual int included_size() const = 0;
        virtual int excluded_size() const = 0;
        virtual int maximum_size() const = 0;

        virtual ~indexer() { }
};

}

#endif
