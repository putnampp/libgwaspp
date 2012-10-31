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
#include "util/lookup_table/lookup_table.h"

namespace util {

//template< class T, class I>
//I LookupTable< T, I >::getIndex( T &t ) throw() {
//    IndexedElement ie( t, -1 );
//    if(( eIter = elements.find( ie ) ) != elements.end() ) {
//        return eIter->second;
//    }
//    throw NotFoundException();
//}
//
//template< class T, class I >
//T LookupTable<T, I>::getElement( I i ) throw() {
//    if( i < 0 || i >= ( int ) elements.size() ) {
//        throw InvalidIndexException();
//    }
//
//    return indexes[i]->first;
//}
//
//template< class T, class I >
//I LookupTable<T, I>::addElement( T &t ) throw() {
//    if(( uint )elements.size() >= max_elements || elements.size() >= elements.max_size() ) {
//        throw OutOfBoundsException();
//    }
//
//    IndexedElement ie( t, ( I ) elements.size() );
//    pair< ElementSetIterator , bool> success = elements.insert( ie );
//    if( success.second ) {
//        indexes.push_back( success.first );
//        return ie.second;
//    }
//    return success.first->second;
//}
//
//template< class T, class I >
//void LookupTable<T, I>::reset() {
//    this->elements.clear();
//    this->indexes.clear();
//}

//template< class T, class I >
//LookupTable<T, I>::~LookupTable() {
//    //dtor
//    reset();
//}

}
