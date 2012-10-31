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
#ifndef TILE_H
#define TILE_H

#include <cstring>

#include "memory/tile_iterator.h"

using namespace std;

namespace util {

template < class D >
class Tile {
    public:
        Tile( int R, int C );
        Tile( const Tile<D> &t );

        D &operator()( int row, int col );

        TileIterator<D> begin_column( int c ) { return TileIterator<D>( mem + c, COL_COUNT ); }
        TileIterator<D> end_column( int c ) { return TileIterator<D>( mem + ELEMENT_COUNT + c ); }

        TileIterator<D> begin_row( int r ) { return TileIterator<D>( mem + r * COL_COUNT ); }
        TileIterator<D> end_row( int r ) { return TileIterator<D>( mem + COL_COUNT ); }

        virtual ~Tile();
    protected:
    private:
        const int ROW_COUNT, COL_COUNT, ELEMENT_COUNT;
        D *mem;
};

template < class D >
Tile<D>::Tile( int R, int C ) : ROW_COUNT( R ), COL_COUNT( C ), ELEMENT_COUNT( R *C ) {
    mem = new D[ ELEMENT_COUNT ];
    memset( mem, 0, ELEMENT_COUNT * sizeof( D ) );
}

template < class D >
Tile<D>::Tile( const Tile<D> &t ) : ROW_COUNT( t.ROW_COUNT ), COL_COUNT( t.COL_COUNT ), ELEMENT_COUNT( t.ELEMENT_COUNT ) {
    mem = new D[ ELEMENT_COUNT ];

    memcpy( mem, t.mem, ELEMENT_COUNT * sizeof( D ) );
}

template < class D >
D &Tile<D>::operator()( int row, int col ) {
    return mem + row * COL_ELE_OFFSET + col;
}

template < class D >
Tile<D>::~Tile() {
    delete [] mem;
}

}

#endif // TILE_H
