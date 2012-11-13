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
#include "genetics/genotype/common_genotype_func.h"

namespace libgwaspp {
namespace genetics {

void printContingencyTable( const contingency_table & ct, ostream & out, bool include_header ) {
    if( include_header ) {
        out << "\tAA\tAa\taa\nBB\t";
    }
    out << dec << ct.n0 << "\t" << ct.n1 << "\t" << ct.n2 << endl;
    if( include_header )
        out << "Bb\t";
    out << ct.n3 << "\t" << ct.n4 << "\t" << ct.n5 << endl;
    if( include_header )
        out << "bb\t";
    out << ct.n6 << "\t" << ct.n7 << "\t" << ct.n8 << endl;
}

void printFrequencyDistribution( const frequency_table & ft, ostream & out, bool include_header ) {
    if( include_header ) {
        out << "AA\tAa\taa\n";
    }
    out << dec << ft.aa << "\t" << ft.ab << "\t" << ft.bb << endl;
}


// Modified version of ones32 found on
// http://aggregate.org/MAGIC/#Population Count (Ones Count)
uint ones16( register uint x ) {
    x -= (( x >> 1 ) & 0x5555 );
    x = ((( x >> 2 ) & 0x3333 ) + ( x & 0x3333 ) );
    x = ((( x >> 4 ) + x ) & 0x0f0f );
    x += ( x >> 8 );
    return ( x & 0x003f );
}

void init_bit_count( byte * _count ) {
    byte * tmp = _count;
    for( uint i = 0; i < 0x10000; ++i, ++tmp ) {
        *tmp = ones16( i );
    }
}

}
}
