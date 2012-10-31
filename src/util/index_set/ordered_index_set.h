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
#ifndef ORDEREDINDEXSET_H
#define ORDEREDINDEXSET_H

#include <string>
#include <vector>
#include <boost/bimap.hpp>

#include "util/index_set/indexer.h"
#include "util/index_set/indexable.h"

using namespace std;

namespace util {

class OrderedIndexSet : public indexer {
    public:
        OrderedIndexSet( indexable * idx, vector< int > & indices );
        OrderedIndexSet( indexable * idx, vector< string > & ids );

        IndexIterator *included_begin() { return new IncludedIndexIterator( current_set, max_elem_count ); }
        IndexIterator *included_end() { return new IncludedIndexIterator( current_set + set_byte_size, -1 ); }

        IndexIterator *excluded_begin() { return new ExcludedIndexIterator( current_set, max_elem_count ); }
        IndexIterator *excluded_end() { return new ExcludedIndexIterator( current_set + set_byte_size, -1 ); }

        virtual ~OrderedIndexSet();
    protected:
        inline void include( int idx, int val );
        void include( const std::vector<int>& idx );
        void include( int *indices, int count );

        void include( const std::string &id );
        void include( const std::vector<std::string>& id );

        inline void exclude( int idx );
        void exclude( const std::vector<int>& idx );
        void exclude( int *indices, int count );

        void exclude( const std::string &id );
        void exclude( const std::vector<std::string>& id );

        void clear();

        indexable *super;
        int *current_set;

        bimap< int, int > order;

        int max_elem_count;
        int included_count, excluded_count;
        int set_byte_size;
};

}

#endif // ORDEREDINDEXSET_H
