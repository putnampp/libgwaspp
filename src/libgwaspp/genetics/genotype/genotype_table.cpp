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
#include "genetics/genotype/genotype_table.h"

namespace libgwaspp {
namespace genetics {

GenotypeTable *GenotypeTable::instance = NULL;

GenotypeTable::GenotypeTable() {
    //ctor
}

GenotypeTable *GenotypeTable::getInstance() {
    if( instance == NULL ) {
        instance = new GenotypeTable();
    }
    return instance;
}

void GenotypeTable::loadGenotypeTable( string &gt_file ) {
}

byte GenotypeTable::findEncodedGenotype( Genotype &gt ) throw() {
    EncodedGenotype eg( gt, -1 );

    if( (gtIterator = genos.find( eg )) != genos.end()) {
        return gtIterator->second;
    }

    throw GenotypeMissingException();
}

Genotype GenotypeTable::getGenotype( byte encoding ) throw() {
    if( encoding < 0 || encoding >= (int) genos.size()) {
        throw libgwaspp::util::InvalidIndexException();
    }

    return index[encoding]->first;
}

byte GenotypeTable::addGenotype( EncodedGenotype &eg ) throw() {
    if( (int)genos.size() >= 255) {
        throw GenotypeTableFullException();
    }
    eg.second = ( byte )genos.size();
    pair<GenotypeSet::iterator, bool> success = genos.insert( eg );
    if( success.second ) {
        index.push_back( success.first );
        return eg.second;
    }
    return success.first->second;
}

void GenotypeTable::reset() {
    index.clear();
    genos.clear();
}

GenotypeTable::~GenotypeTable() {
    //dtor
    reset();
}

}
}
