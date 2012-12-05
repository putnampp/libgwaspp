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
#ifndef COMMON_GENOTYPE_FUNC_H
#define COMMON_GENOTYPE_FUNC_H

#include <iostream>
#include <fstream>
#include <cassert>

#define MATHLIB_STANDALONE
#include "Rmath.h"

#include "genetics/genotype/common_genotype.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

void printFrequencyDistribution( const frequency_table & ft, ostream & out, bool include_header = false );
void printContingencyTable( const contingency_table & ct, ostream & out, bool include_header = false );

ostream &operator<<(ostream &out, contingency_table & ct );
ostream &operator<<(ostream &out, frequency_table &ft );
ostream &operator<<(ostream &out, marginal_information &mi );

uint ones16( register uint i );
int init_bit_count( byte * bit_count );

inline int PopCount( ushort v ) {
    static byte bit_count16[ 0x10000 ];
    static int init = init_bit_count( bit_count16 );
    return bit_count16[ v ];
}

/*
inline int PopCount( ulong v ) {
    static int count;
    count =  PopCount( (ushort) v );
    v >>= 16;
    count += PopCount( (ushort) v );
    v >>= 16;
    count += PopCount( (ushort) v );
    v >>= 16;
    count += PopCount( (ushort) v );
    return count;
}

inline int PopCount( uint v ) {
    static int count;
    count = PopCount( (ushort) v );
    v >>= 16;
    count += PopCount( (ushort) v );
    return count;
}
*/

inline int PopCount( register uint v ) {
    return PopCount( (ushort) v ) + PopCount( (ushort) (v >> 16) );
}

inline int PopCount( register ulong v ) {
    return PopCount( (ushort) v ) + PopCount( (ushort) (v >> 16) ) + PopCount( (ushort) (v >> 32) ) + PopCount( (ushort) (v >> 48) );
}

#ifndef IncrementFrequencyValueStream
#define IncrementFrequencyValueStream( ft, l_aa, l_ab, l_bb )       \
                            ft.aa += PopCount( l_aa );              \
                            ft.ab += PopCount( l_ab );              \
                            ft.bb += PopCount( l_bb );
#endif

#ifndef AddToContingencyStream
#define AddToContingencyStream( n, a, b )    \
    n += PopCount( (a & b ) );
#endif

void computeMarginalInformation( const frequency_table & _cases, const frequency_table & _ctrls, uint nIndivids, marginal_information & m);

}
}

#endif
