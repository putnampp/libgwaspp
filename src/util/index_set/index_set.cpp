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
#include "util/index_set/index_set.h"

namespace util {

//IndexedSet::IndexedSet( indexable *idx, bool all_in ) : super( idx ),  max_elem_count( idx->size() ), included_count( 0 ), excluded_count( 0 ) {
//    set_byte_size = byte_count( max_elem_count );
//
//    // Allocate 1-bit for each indexable element from the indexable super set
//    current_set = new byte[ set_byte_size + 1];
//    current_set[set_byte_size] = 0;
//
//    if( all_in ) {
//        // Assume that every element is initially included
//        includeAll();
//    } else {
//        // Assume that every element is initially excluded
//        excludeAll();
//    }
//}
//
//IndexedSet::IndexedSet( const IndexedSet & idxset) :
//        super( idxset.super ), max_elem_count( idxset.max_elem_count ),
//        included_count( idxset.included_count), excluded_count(idxset.excluded_count), set_byte_size( idxset.set_byte_size) {
//    current_set = new byte[ set_byte_size + 1];
//    current_set[set_byte_size] = 0;
//
//    memcpy( current_set, idxset.current_set, set_byte_size);
//}

IndexedSet::IndexedSet( indexable *idx, const std::vector<int> &indices ) : super( idx ),  max_elem_count( idx->size() ), included_count( 0 ), excluded_count( max_elem_count ) {
    set_byte_size = byte_count( max_elem_count );

    this->current_set = new byte[ set_byte_size + 1 ];
    this->current_set[ set_byte_size ] = 0;
    memset( this->current_set, ( char ) 0, set_byte_size );


    this->order = new int[( int )indices.size() + 1];
    this->order[( int )indices.size()] = 0;
    memset( order, ( char ) - 1, ( int )indices.size() * sizeof( int ) );

    include( indices );
}

IndexedSet::IndexedSet( indexable *idx, const std::vector<std::string> & ids ) : super( idx ),  max_elem_count( idx->size() ), included_count( 0 ), excluded_count( max_elem_count ) {
    set_byte_size = byte_count( max_elem_count );

    this->current_set = new byte[ set_byte_size + 1 ];
    this->current_set[ set_byte_size ] = 0;
    memset( this->current_set, ( char ) 0, set_byte_size );

    this->order = new int[( int )ids.size() + 1];
    this->order[( int )ids.size()] = 0;
    memset( order, ( char ) - 1, ( int )ids.size() * sizeof( int ) );

    include( ids );
}

int IndexedSet::orderOf( const std::string &id ) {
    int idx = ( *super )( id );
    int ord = -1;
    int *order_ptr = order;

    // linear search is bad!!
    for( int i = 0; i < included_count && ord < 0; ++i, ++order_ptr ) {
        ord = ( *order_ptr == idx ) ? i : -1;
    }

    return ord;
}

int IndexedSet::indexOf( int _ord ) {
    return order[_ord];
}

std::string IndexedSet::getIDAtOrderedIndex( int _ord ) {
    assert( order[ _ord ] > -1 && order[_ord] < max_elem_count );
    return ( *super )( order[ _ord ] );
}

void IndexedSet::includeAll() {
    memset( current_set, ( char ) - 1, set_byte_size );
    included_count = max_elem_count;
    excluded_count = 0;
}

void IndexedSet::include( int idx, int ord ) {
    if( !bitset_isset( current_set, idx ) ) {
        bitset_set_idx( current_set, idx );
        ++included_count;
        --excluded_count;

        order[ord] = idx;
    }
}

void IndexedSet::include( int *indices, int len ) {
    int *tmp = indices;
    for( int i = 0; i < len; ++i, ++tmp ) {
        include( *tmp, i );
    }
}

void IndexedSet::include( const std::string &id, int at ) {
    int idx = ( *super )( id );
    include( idx, at );
}

void IndexedSet::include( const std::vector<int>& indices ) {
    int ord = 0;
    for( std::vector<int>::const_iterator it = indices.begin(), end = indices.end(); it != end; ++it, ++ord ) {
        include( *it, ord );
    }
}

void IndexedSet::include( const std::vector<std::string> & ids ) {
    int idx = 0, at = 0;
    for( std::vector<std::string>::const_iterator it = ids.begin(); it != ids.end(); ++it, ++at ) {
        idx = ( *super )( *it );
        include( idx, at );
    }
}

void IndexedSet::exclude( int at ) {
    if( bitset_isset( current_set, order[at] ) ) {
        bitset_clear_idx( current_set, at );
        --included_count;
        ++excluded_count;

        order[ at ] = -1;
    }
}

void IndexedSet::exclude( int *indices, int len ) {
    int *tmp = indices;
    for( int i = 0; i < len; ++i, ++tmp ) {
        exclude( *tmp );
    }
}

void IndexedSet::exclude( const std::string &id ) {
    int idx = ( *super )( id );
    exclude( idx );
}

void IndexedSet::exclude( const std::vector<int>& indices ) {
    if(( int )indices.size() != super->size() ) {
        for( std::vector<int>::const_iterator it = indices.begin(); it != indices.end(); ++it ) {
            exclude( *it );
        }
    } else {
        excludeAll();
    }
}

void IndexedSet::exclude( const std::vector<std::string>& ids ) {
    if(( int ) ids.size() != super->size() ) {
        int idx;
        for( std::vector<std::string>::const_iterator it = ids.begin(); it != ids.end(); ++it ) {
            idx = ( *super )( *it );
            exclude( idx );
        }
    } else {
        excludeAll();
    }
}

void IndexedSet::excludeAll() {
    // Assume that every element is initially included
    memset( current_set, ( char )0, set_byte_size );
    included_count = 0;
    excluded_count = max_elem_count;
}

IndexedSet::~IndexedSet() {
    delete [] current_set;
    delete [] order;
}

}
