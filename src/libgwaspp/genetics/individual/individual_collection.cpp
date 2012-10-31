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
#include "genetics/individual/individual_collection.h"

namespace libgwaspp {
namespace genetics {

Individual *IndividualCollection::findOrCreateIndividual( string &id ) {
    IndLookupSet::iterator ptr_iter;

    if(( ptr_iter = individ_lookup.find( id ) ) == individ_lookup.end() ) {
        Individual *ind = new Individual( id, ptree );

        individ_lookup[ id ] = (int) individs.size();
        individs.push_back( ind );

        return ind;
    } else {
        return individs[ptr_iter->second];
    }
}

Individual *IndividualCollection::getIndividualAt( int idx ) {
    return individs[ idx ];
}


int IndividualCollection::operator()( const string &id ) {
    IndLookupSet::const_iterator ptr_iter;
    if(( ptr_iter = individ_lookup.find( id ) ) == individ_lookup.end() ) {
        Individual *ind = new Individual( id, ptree );

        individ_lookup[ id ] = (int) individs.size();
        individs.push_back( ind );

        return (int) individs.size() - 1;
    }

    return ptr_iter->second;
}

string IndividualCollection::operator()( int idx ) const {
    return individs[ idx ]->getID();
}


void IndividualCollection::reset() {
    for( vector< Individual * >::iterator it = individs.begin(); it != individs.end(); it++ ) {
        if( *it != NULL ) {
            delete *it;
        }
    }
    individ_lookup.clear();
    individs.clear();
    ptree->reset();
}

IndividualCollection::~IndividualCollection() {
    //dtor

    reset();

    delete ptree;
}

}
}
