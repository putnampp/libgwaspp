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
#ifndef TABLE_H
#define TABLE_H

#include <iostream>
#include <vector>

#include "util/index_set/index_set.h"

using namespace std;

namespace util {

/**
        A table is simply a way of representing ordered data in as
        contiguous sets of contiguous data elements. From an abstract
        perspective, a table only needs to know the dimensions
        for defining the contiguous regions.
*/

template < class D >
class Table {
    public:
        Table( indexer *r, indexer *c ) : rows( r ), columns( c ), max_row( r->included_size() ), max_column( c->included_size() ), data( NULL ) { }

        virtual D operator()( int r, int c ) { return data[ r * max_column + c ]; }

        int row_size() const { return max_row; }
        int column_size() const { return max_column; }

        string rowID( int row_idx ) const { return rows->getIDAtOrderedIndex( row_idx ); }
        string columnID( int column_idx ) const { return columns->getIDAtOrderedIndex( column_idx ); }

        const indexer * getRowSet() { return rows; }
        const indexer * getColumnSet() { return columns; }

        virtual ~Table() {
            delete [] data;
        }
    protected:
        Table( int maxR, int maxC ) : rows( NULL ), columns( NULL ), max_row( maxR ), max_column( maxC ), data( NULL ) { }
        virtual void initialize() = 0;

        indexer *rows,  *columns;
        int max_row, max_column;
        D *data;
};

}

#endif // TABLE_H
