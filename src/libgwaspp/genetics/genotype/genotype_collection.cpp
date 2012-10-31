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
#include "genetics/genotype/genotype_collection.h"

namespace libgwaspp {
namespace genetics {

GenotypeCollection::GenotypeCollection() {
    //ctor
}

const Genotype *GenotypeCollection::findOrCreateGenotype( string &gt ) {
    LookupTable::iterator it = lookup.find( &gt );

    if( it == lookup.end() ) {
        // gt does not exist
        Genotype *g = new SNPGenotype( gt );
        genotypes.push_back( g );
        lookup.insert( pair< const string *, GenotypeID>( &g->geno, ( GenotypeID )genotypes.size() - 1 ) );
        return g;
    }

    return genotypes[ it->second ];
}

GenotypeID GenotypeCollection::findOrCreateGenotypeID( string &gt ) {
    if( last_gt == gt ) {
        return last_gt_id;
    }

//    LookupTable::iterator it = lookup.find( &gt );
//
//    if( it == lookup.end() ) {
//        Genotype *g = new SNPGenotype(gt);
//        genotypes.push_back( g );
//        pair< LookupTable::iterator, bool > success =
//            lookup.insert( pair< const string *, GenotypeID>( &g->geno, (GenotypeID) genotypes.size() - 1) );
//
//        last_gt = *(success.first->first);
//        last_gt_id = success.first->second;
//
//        return last_gt_id;
//    }
//
//    last_gt = *(it->first);
//    last_gt_id = it->second;

    for( LookupTable::iterator it = lookup.begin(); it != lookup.end(); it++ ) {
        if( *it->first == gt ) {
            last_gt = gt;
            last_gt_id = it->second;
            return last_gt_id;
        }
    }

    Genotype *g = new SNPGenotype( gt );
    genotypes.push_back( g );
    pair< LookupTable::iterator, bool > success =
        lookup.insert( pair< const string *, GenotypeID>( &g->geno, ( GenotypeID ) genotypes.size() - 1 ) );

    last_gt = *( success.first->first );
    last_gt_id = success.first->second;

    return last_gt_id;
}

GenotypeID GenotypeCollection::getGenotypeID( string &gt ) {
    LookupTable::iterator it = lookup.find( &gt );

    if( it == lookup.end() ) {
        return -1;
    }
    return it->second;
}

GenotypeCollection::~GenotypeCollection() {
    //dtor

    for( vector< Genotype * >::iterator it = genotypes.begin(); it != genotypes.end(); it++ ) {
        delete *it;
    }

    lookup.clear();
    genotypes.clear();
}

}
}
