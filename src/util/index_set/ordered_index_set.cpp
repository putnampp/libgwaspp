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
#include "util/index_set/ordered_index_set.h"

namespace util {

OrderedIndexSet::OrderedIndexSet( indexable *idx, const vector< int > & indices ) : super( idx ), max_elem_count( idx->size() ), included_count( 0 ), excluded_count( 0 )  {
    //ctor
    current_set = new int[ max_elem_count + 1];
    current_set[ max_elem_count ] = 0;

    clear();
    include( indices );
}

void OrderedIndexSet::include( int idx, int at ) {
    assert( current_set[ idx ] == -1 );

    current_set[ idx ] = at;
    ++included_count;
    --excluded_count;
}

void OrderedIndexSet::include( const std::vector<int>& idx ) {
    int at = 0;
    for( vector<int>::const_iterator it = idx.begin(); it != idx.end(); ++it, ++at ) {
        include( *it, at );
    }
}

void OrderedIndexSet::include( int *indices, int count ) {
    for( int i = 0; i < count; ++i ) {
        include( indices[i], i );
    }
}

void OrderedIndexSet::include( const std::string &id, int at ) {
    include(( *super )( id ), at );
}

void OrderedIndexSet::include( const std::vector<std::string>& id ) {
    int at = 0;
    for( vector<string>::const_iterator it = id.begin(); it != id.end(); ++it, ++at ) {
        include(( *super )( *it ), at );
    }
}

void OrderedIndexSet::exclude( int idx ) {
    assert( current_set[ idx ] >= 0 );

    current_set[ idx ] = -1;
    ++excluded_count;
    --included_count;
}

void OrderedIndexSet::exclude( const std::vector<int>& idx ) {
    for(vector<int>::const_iterator it = idx.begin(); it != idx.end(); ++it) {
        exclude( *it );
    }
}

void OrderedIndexSet::exclude( int *indices, int count ) {
    for( int i = 0; i < count; ++count ) {
        exclude( indices[i], i );
    }
}

void OrderedIndexSet::exclude( const std::string &id ) {
    exclude( (*super)(id) );
}

void OrderedIndexSet::exclude( const std::vector<std::string>& id ) {
    for( vector<string>::const_iterator it = id.begin(); it != id.end(); ++it) {
        exclude( (*super)(*it) );
    }
}

OrderedIndexSet::clear() {
    memset( current_set, ( char ) - 1, max_elem_count * sizeof( int ) );
    excluded_count = max_elem_count;
    included_count = 0;
}

OrderedIndexSet::~OrderedIndexSet() {
    //dtor
    delete [] current_set;
}

}
