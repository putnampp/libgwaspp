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
#include "marker_collection.h"

namespace libgwaspp {
namespace genetics {

MarkerCollection::MarkerCollection() {
    //ctor
    alleles = new AlleleCollection();
    chromosomes = new ChromosomeCollection();
}

const Marker *MarkerCollection::findMarker( string &id ) const {
    LookupMap::const_iterator it = marker_lookup.find( &id );

    if( it == marker_lookup.end() ) {
        return NULL;
    }
    return markers[ it->second ];
}

int MarkerCollection::operator()( const string &id ) {
    LookupMap::const_iterator it = marker_lookup.find( &id );

    if( it == marker_lookup.end() ) {
        return -1;
    }
    return it->second;
}

string MarkerCollection::operator()( int idx ) const {
    return markers[ idx ]->getID();
}

const Marker *MarkerCollection::createMarker( string &id, string &chrom, uint start, uint end, double gPos, string &alls ) {
    const Marker *m = findMarker( id );

    if( m == NULL ) {
        // create new Marker
        Marker *m2 = new Marker( id, chromosomes->createChromosome( chrom ), start, end, gPos, alleles->findOrCreateAllele(alls) );

        markers.push_back( m2 );
        marker_lookup.insert( pair< const string *, int>( m2->getIDPtr(), (int) markers.size() - 1));

        return m2;
    }

    return m;
}

MarkerCollection::~MarkerCollection() {
    //dtor
    for( vector< Marker * >::iterator it = markers.begin(); it != markers.end(); it++ ) {
        delete *it;
    }

    markers.clear();
    marker_lookup.clear();

    delete alleles;
    delete chromosomes;
}

}
}
