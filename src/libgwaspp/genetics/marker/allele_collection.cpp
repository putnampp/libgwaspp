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
#include "genetics/marker/allele_collection.h"

namespace libgwaspp {
namespace genetics {

const AlleleForm * AlleleCollection::findOrCreateAllele( string & alls, char delim ) {
    ulong hash = AlleleForm::computeHashID( alls );

    map< ulong, ALLELE_INDEX >::iterator it = key_lookup.find( hash );

    if( it == key_lookup.end() ) {
        // new Allele

        AlleleForm *a = new AlleleForm( alls, delim );

        pair< map< ulong, ALLELE_INDEX>::iterator, bool> success =
            key_lookup.insert( pair< ulong, ALLELE_INDEX >( a->getID(), (ALLELE_INDEX) alleles.size()) );
        alleles.push_back( a );

        return a;
    }

    return alleles[ it->second ];
}

ALLELE_INDEX AlleleCollection::findOrCreateAlleleIndex( string & alls, char delim ) {
    ulong hash = AlleleForm::computeHashID( alls );

    map< ulong, ALLELE_INDEX >::iterator it = key_lookup.find( hash );

    if( it == key_lookup.end() ) {
        // new Allele
        AlleleForm *a = new AlleleForm( alls, delim );

        pair< map< ulong, ALLELE_INDEX >::iterator, bool> success =
            key_lookup.insert( pair< ulong, ALLELE_INDEX >( a->getID(), (ALLELE_INDEX) alleles.size() ) );
        alleles.push_back( a );

        return success.first->second;
    }

    return it->second;
}

ALLELE_INDEX AlleleCollection::findIndexOf( string &alls ) {
    ulong hash = AlleleForm::computeHashID( alls );
    map< ulong, ALLELE_INDEX >::iterator it = key_lookup.find( hash );

    if( it == key_lookup.end() ) {
        return -1;
    }

    return it->second;
}

AlleleCollection::~AlleleCollection() {
    //dtor
    for( vector< AlleleForm * >::iterator it = alleles.begin(); it != alleles.end(); it++ ) {
        delete * it;
    }

    alleles.clear();
}

}
}
