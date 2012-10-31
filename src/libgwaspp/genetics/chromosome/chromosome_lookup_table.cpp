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
#include "genetics/chromosome/chromosome_lookup_table.h"


namespace libgwaspp {
namespace genetics {

ChromosomeLookupTable *ChromosomeLookupTable::instance = NULL;

ChromosomeLookupTable::ChromosomeLookupTable() : LookupTable<Chromosome>() {
    //ctor
}

ChromosomeLookupTable *ChromosomeLookupTable::getInstance() {
    if( instance == NULL ) {
        instance = new ChromosomeLookupTable();
    }
    return instance;
}

byte ChromosomeLookupTable::getIndex( string &name) throw() {
    if( (aliasIter = alias.find(name)) != alias.end() ) {
        return aliasIter->second->second;
    }
    throw NotFoundException();
}

byte ChromosomeLookupTable::addElement( Chromosome &chrom) throw() {
    if( (uint) elements.size() >= this->max_elements || elements.size() >= elements.max_size()) {
        throw OutOfBoundsException();
    }
    IndexedElement ie( chrom, (byte) elements.size() );
    pair< ElementSetIterator, bool> success = elements.insert( ie );
    if( success.second ) {
        alias[ chrom.getName() ] = success.first;
        return ie.second;
    }
    return success.first->second;
}

byte ChromosomeLookupTable::addAlias( Chromosome &c, string &altName ) throw() {
    if( (aliasIter = alias.find( c.getName() ) ) != alias.end() ) {
        alias[ altName ] = aliasIter->second;
        return aliasIter->second->second;
    }
    throw NotFoundException();
}

byte ChromosomeLookupTable::addAlias( byte idx, string &altName ) throw() {
    if( idx >= (byte) elements.size()) {
        throw InvalidIndexException();
    }

    alias[ altName ] = indexes[ idx ];

    return idx;
}

void ChromosomeLookupTable::reset() {
    LookupTable<Chromosome>::reset();

    alias.clear();
}

ChromosomeLookupTable::~ChromosomeLookupTable() {
    //dtor
}

}
}
