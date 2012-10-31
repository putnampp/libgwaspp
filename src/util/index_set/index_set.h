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
#ifndef INDEXSET_H
#define INDEXSET_H

#include <iostream>
#include <cstring>
#include <vector>
#include <iterator>
#include <cassert>

#include "common.h"

#include "util/index_set/index_iterator.h"
#include "util/index_set/indexable.h"
#include "util/index_set/indexer.h"

using namespace std;

namespace util {


class ExcludedIndexIterator : virtual public IndexIterator {
    public:
        ExcludedIndexIterator( byte *p, int max_idx ) : ptr( p ), val( *p ), current_idx( -1 ), maximum( max_idx ) {
            int offset;
            while( current_idx < maximum ) {
                if( !( offset = bitset_bit_offset( ++current_idx ) ) ) { val = ( byte ) * ptr; }
                if( !is_bit_set( val, offset ) )    break;
            }
        }

        ExcludedIndexIterator( const ExcludedIndexIterator &eii ) : ptr( eii.ptr ), val( eii.val ), current_idx( eii.current_idx ), maximum( eii.maximum ) {}

        IndexIterator *operator++() {
            int offset;
            while( current_idx < maximum ) {
                if( !( offset = bitset_bit_offset( ++current_idx ) ) ) { val = ( byte ) * ptr; }
                if( !is_bit_set( val, offset ) )   break;
            }
            return this;
        }

        inline bool operator==( const ExcludedIndexIterator &rhs ) const { return ptr == rhs.ptr; }
        inline bool operator!=( const ExcludedIndexIterator &rhs ) const { return !operator==( rhs ); }

        int operator*() { return current_idx; }

        bool is_done() { return current_idx >= maximum; }

    protected:
        byte *ptr;
        byte val;
        int current_idx, maximum;
};

class IncludedIndexIterator : virtual public IndexIterator {
    public:
        IncludedIndexIterator( byte *p, int max_idx ) : ptr( p ), val( *p ), current_idx( -1 ), maximum( max_idx ) {
            int offset;
            while( current_idx < maximum ) {
                if( !( offset = bitset_bit_offset( ++current_idx ) ) ) { val = *ptr; }
                if( is_bit_set( val, offset ) )    break;
            }
        }

        IndexIterator *operator++() {
            int offset;
            while( current_idx < maximum ) {
                if( !( offset = bitset_bit_offset( ++current_idx ) ) ) { val = *ptr; }
                if( is_bit_set( val, offset ) )    break;
            }
            return this;
        }

        inline bool operator==( const IncludedIndexIterator &rhs ) const { return ptr == rhs.ptr; }
        inline bool operator!=( const IncludedIndexIterator &rhs ) const { return !operator==( rhs ); }

        int operator*() { return current_idx; }

        bool is_done() { return current_idx >= maximum; }
    protected:
        byte *ptr;
        byte val;
        int current_idx, maximum;
};

class IndexedSet : public indexer {
    public:
//        IndexedSet( indexable *sup, bool all_in = true );
        IndexedSet( indexable *sup, const std::vector<int>& indices );
        IndexedSet( indexable *sup, const std::vector<std::string> & ids );

        bool operator==( const IndexedSet &rhs ) { return this->super == rhs.super; }

        int orderOf( const std::string &id );
        int indexOf( int _ord );

        std::string getIDAtOrderedIndex( int _ord );

        IndexIterator *included_begin() { return new IncludedIndexIterator( current_set, max_elem_count ); }
        IndexIterator *included_end() { return new IncludedIndexIterator( current_set + set_byte_size, -1 ); }

        IndexIterator *excluded_begin() { return new ExcludedIndexIterator( current_set, max_elem_count ); }
        IndexIterator *excluded_end() { return new ExcludedIndexIterator( current_set + set_byte_size, -1 ); }

        int maximum_size() const { return max_elem_count; }
        int included_size() const { return included_count; }
        int excluded_size() const { return excluded_count; }

        virtual ~IndexedSet();

    protected:
        inline void include( int idx, int at );
        void include( const std::vector<int>& idx );
        void include( int *indices, int count );

        void include( const std::string &id, int at );
        void include( const std::vector<std::string>& id );

        void includeAll();

        inline void exclude( int idx );
        void exclude( const std::vector<int>& idx );
        void exclude( int *indices, int count );

        void exclude( const std::string &id );
        void exclude( const std::vector<std::string>& id );

        void excludeAll();

        indexable * super;
        byte *current_set;

        int *order;

        int max_elem_count;
        int included_count, excluded_count;
        int set_byte_size;
};

}

#endif // INDEXSET_H
