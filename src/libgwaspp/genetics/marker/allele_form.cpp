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
#include "genetics/marker/allele_form.h"

namespace libgwaspp {
namespace genetics {

AlleleForm::AlleleForm( string &alleles, char delim ) {
    //ctor
    initialize( alleles, delim );
}

ulong AlleleForm::computeHashID( string &id ) {
    ulong hash = 0;

    for( string::iterator it = id.begin(); it != id.end(); it++ ) {
        hash = ( hash << 5 ) ^( hash >> 22 ) ^ *it;
    }

    return hash;
}

void AlleleForm::initialize( string &alls, char delim ) {
    this->id = computeHashID( alls );
    istringstream iss( alls );
    string tok;

    if( delim == 0 ) {
        char val;
        while( !iss.eof() ) {
            val = ( char )iss.get();
            alleles.push_back( string( 1, val ) );
            iss.peek();
        }
    } else {
        while( getline( iss, tok, delim ) ) {
            alleles.push_back( tok );
        }
    }
}

void AlleleForm::addAllele( string &a ) {
    alleles.push_back( a );
}

bool AlleleForm::operator==( const AlleleForm &af ) {
    if( this->alleles.size() == af.alleles.size() ) {
        for( vector< string >::const_iterator it = this->alleles.begin(), it2 = af.alleles.begin(); it != this->alleles.end(); it++, it2++ ) {
            if( *it != *it2 ) return false;
        }
        return true;
    }
    return false;
}



AlleleForm::~AlleleForm() {
    //dtor
    alleles.clear();
}

}
}
