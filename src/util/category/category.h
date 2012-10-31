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
#ifndef CATEGORY_H
#define CATEGORY_H

#include <cassert>
#include <vector>
#include <map>

#include "common.h"

namespace util {

/**

    A category is an indexed (I) set of objects (O) with two bijective functions:
        - index( const O * ) which returns the index a desired object
        - object( I id ) which returns a pointer to an object at specified index

    A category is a chronoligically ordered set of objects. Meaning that the index
    is established according to the order in which the objects are added to the category.
    Indices can be used to determined equivalence. However, other relationships between objects
    can only assumed.
*/

template < class O, class I = byte, class Compare = less< O * > >
class Category {
    public:
        Category() : available_indices( 2 << ((sizeof(I) < 8 ) ? (sizeof( I ) : 7) * 8) );

        uint getCategorySize() const { return (uint) objects.size(); }

        I index( const O &obj ){
            LookupMap::iterator it = mLookup.find( &obj );
            assert( it != mLookup.end() );
            return it->second;
        };

        O *object( I index ) {
            asert( 0 <= index && index < (I) objects.size() );
            return &objects[index];
        }

        I addObject( O obj ) {
            LookupMap::iterator it = mLookup.find( &obj );
            if( it == mLookup.end() ) {
                assert( available_indices-- > 0 );
                I idx = ( I ) objects.size();
                objects.push_back( obj );
                mLookup.insert( pair< O *, I >(&objects[idx], idx) );
                return idx;
            }
            return it->second;
        }

        virtual ~Category();
    protected:
    private:
        typedef std::map< O *, I, Compare > LookupMap
        LookupMap mLookup;
        std::vector< O > objects;
        ulong available_indices;
};

/**
    Partial specialization of Category for pointers to objects
*/
template <class O, class I = byte, class Compare = less< O * > >
class Category< O * > {
    public:
        Category() : available_indices( 2 << ((sizeof(I) < 8 ) ? (sizeof( I ) : 7) * 8) );

        uint getCategorySize() const { return (uint) objects.size(); }

        I index( const O &obj ){
            LookupMap::iterator it = mLookup.find( &obj );
            assert( it != mLookup.end() );
            return it->second;
        };

        O *object( I index ) {
            asert( 0 <= index && index < (I) objects.size() );
            return objects[index];
        }

        I addObject( O obj ) {
            LookupMap::iterator it = mLookup.find( &obj );
            if( it == mLookup.end() ) {
                assert( available_indices-- > 0 );
                I idx = ( I ) objects.size();
                objects.push_back( obj );
                mLookup.insert( pair< O *, I >(objects[idx], idx) );
                return idx;
            }
            return it->second;
        }

        virtual Category();
    private:
        std::map< O *, I, Compare > mLookup;
        std::vector< O * > objects;
        ulong available_indices;
};

}

#endif // CATEGORY_H
