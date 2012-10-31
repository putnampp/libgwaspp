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
#include "genetics/chromosome/chromosome_table.h"

namespace libgwaspp {
namespace genetics {

ChromosomeTable *ChromosomeTable::instance = NULL;

ChromosomeTable::ChromosomeTable() {
    //ctor
}

ChromosomeTable *ChromosomeTable::getInstance() {
    if( instance == NULL ) {
        instance = new ChromosomeTable();
    }

    return instance;
}

Chromosome ChromosomeTable::getChromosome( byte code ) throw() {
    if( code < 0 || code > (int) index.size()) {
        throw util::InvalidIndexException();
    }
    return index[ code ]->first;
}

byte ChromosomeTable::findChromosomeEncoding( string &c ) throw() {
    if( (aliasIter = chrom_aliases.find( c ) ) != chrom_aliases.end() ) {
        return aliasIter->second->second;
    }

    throw ChromosomeNotFoundException();
}

byte ChromosomeTable::findChromosomeEncoding( const Chromosome &chrom ) throw() {
    EncodedChromosome ec( chrom, ( byte )255 );
    if(( chromIter = chroms.find( ec ) ) != chroms.end() ) {
        return chromIter->second;
    }

    throw ChromosomeNotFoundException();
}



void ChromosomeTable::loadChromosomeTable( string &chromTable ) {
    reset();
}

byte ChromosomeTable::addChromosome( Chromosome &c ) throw() {
    if( chroms.size() >= 255 ) {
        throw ChromosomeTableFullException();
    }
    EncodedChromosome ec( c, ( byte ) chroms.size() );

    pair<ChromosomeSet::iterator, bool> success = chroms.insert( ec );
    if( success.second ) {
        index.push_back( success.first );
        chrom_aliases[ec.first.getName()] = success.first;
        return ec.second;
    }

    return success.first->second;
}

void ChromosomeTable::addChromosomeAlias( Chromosome &c, string &name ) throw() {
    EncodedChromosome ec( c, ( byte )255 );

    if(( chromIter = chroms.find( ec ) ) != chroms.end() ) {
        chrom_aliases[ name ] = chromIter;
    } else {
        throw ChromosomeNotFoundException();
    }
}

void ChromosomeTable::reset() {
    chroms.clear();
    index.clear();
    chrom_aliases.clear();
}

ChromosomeTable::~ChromosomeTable() {
    //dtor

    reset();
}

}
}
