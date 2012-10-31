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
#ifndef MARKERCOLLECTION_H
#define MARKERCOLLECTION_H

#include <iostream>
#include <map>
#include <vector>

#include "libgwaspp.h"
#include "genetics/marker/marker.h"
#include "genetics/marker/allele_collection.h"
#include "genetics/chromosome/chromosome_collection.h"
#include "genetics/chromosome/chromosome.h"
#include "util/index_set/indexable.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

class MarkerCollection : public util::indexable {
    public:
        typedef map< const string *, int, StringPtrComparer > LookupMap;
        MarkerCollection();

        const Marker *findMarker( string &id ) const;
        const Marker *createMarker( string &id, string &chrom, uint start, uint end, double gPos, string &alleles );

        const Marker *getMarkerAt( int idx ) const { return markers[ idx ]; }
        const AlleleForm *getAlleleFormAt( int idx ) const { return alleles->getAlleleAt(idx); }

        int operator()( const string &id );
        string operator()( int idx ) const;

        int size() const { return (int)markers.size(); }

        int getAlleleFormCount() const { return alleles->getAlleleCount(); }
        int getMarkerCount() const { return (int) markers.size(); }

        int getMarkerAlleleCount( string &id );
        int getMarkerAlleleCount( const Marker *m ) const { return m->getAlleles()->getAlleleCount(); }

//        vector< AlleleForm * >::const_iterator alleles_begin() { return alleles.begin(); }
//        vector< AlleleForm * >::const_iterator alleles_end() { return alleles.end(); }

        vector< Marker * >::iterator marker_begin() { return markers.begin(); }
        vector< Marker * >::iterator marker_end() { return markers.end(); }

        virtual ~MarkerCollection();
    protected:
    private:
        LookupMap marker_lookup;
        vector< Marker * > markers;

        AlleleCollection *alleles;
        ChromosomeCollection *chromosomes;

};

}
}

#endif // MARKERCOLLECTION_H
