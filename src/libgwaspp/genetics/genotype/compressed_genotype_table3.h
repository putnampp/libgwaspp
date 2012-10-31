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
#ifndef COMPRESSEDGENOTYPETABLE3_H
#define COMPRESSEDGENOTYPETABLE3_H

#include <iostream>
#include <cmath>

#include "common.h"
#include "util/index_set/indexer.h"
#include "genetics/genotype/geno_table.h"
#include "genetics/genotype/common_genotype.h"

using namespace std;
using namespace util;

namespace libgwaspp {
namespace genetics {

/**
    byte compressed genotype
    2 bits per Genotype
    1 block per Record to identify which genotype code is used as "UNKNOWN" for record

    Little-endian bit ordering
    BLOCK INDEX (ushort):          |        HEADER       |       0       |       1       |  ...
    GENOTYPE INDEX (half-byte):    | 0  |  1  |  2  |  3 |0|1|2|3|4|5|6|7| | | | | | | | |  ...
    CODE INDEX (bit):              |15  |11   |7    |3   |15     |7      |15    8|7     0|  ...

    Genotype Code:
    0000 - AA
    0001 - AC
    0010 - AG
    0011 - AT
    0100 - CA
    0101 - CC
    0110 - CG
    0111 - CT
    1000 - GA
    1001 - GC
    1010 - GG
    1011 - GT
    1100 - TA
    1101 - TC
    1110 - TG
    1111 - TT

    Header:
    0 - Genotype ORDER
    1 - Genotype 1 - (Dom Homo | Het | Rec Homo)
    2 - Genotype 2 - (Dom Homo | Het | Rec Homo)
    3 - Genotype 3 - (Dom Homo | Het | Rec Homo)

    Genotype ORDER (Revised):
    0   ->  { UNSET }
    1   ->  { 1 => Homo; 2 => unset; 3 => unset }
    2   ->  { 1 => unset; 2=> Het; 3 => unset }
    3   ->  { 1 => Homo; 2 => Het; 3 => unset }
    4   ->  Not used;
    5   ->  { 1 => Homo; 2 => unset; 3 => Homo }
    6   ->  Not used;
    7   ->  { 1 => Homo; 2 => Het; 3 => Homo }

    (Not Yet Implemented)
    8   ->  { 1 => Dom Homo; 2 => unset; 3 => unset }
    9   ->  { 1 => Rec Homo; 2 => unset; 3 => unset }
    10  ->  { 1 => Dom Homo; 2 => unset; 3 => Rec Homo }
    11  ->  { 1 => Rec Homo; 2 => unset; 3 => Dom Homo }
    12  ->  { 1 => Dom Homo; 2 => Het; 3 => Rec Homo }
    13  ->  { 1 => Rec Homo; 2 => Het; 3 => Dom Homo }
    14  ->  { 1 => Dom Homo; 2 => Het; 3 => unset }
    15  ->  { 1 => Rec Homo; 2 => Het; 3 => unset }

    Unknown Genotype value = 0xFFFF
*/

//#define HammingWeight( x ) bit_count16[ x & 0x0000FFFF ] + bit_count16[ (x >> 32) & 0x0000FFFF]

#if PROCESSOR_WORD_SIZE == 64

#define U_HALF_MASK 0xAAAAAAAAAAAAAAAA
#define BLOCKS_PER_UNIT 4
#define TAIL_END 31

#ifndef AddToContingency
#define AddToContingency( n, t, a, b )                                      \
    t = (a & b);                                    \
    t = ((t | ((t & 0x5555000055550000) >> 15)));   \
    n += bit_count16[ t & 0x0000FFFF ] + bit_count16[ (t >> 32) & 0x0000FFFF ];
#endif

#ifndef AddToDistribution
#define AddToDistribution(n, a )                                            \
    a = ((a | ((a & 0x5555000055550000) >> 15)))    \
        n += bit_count16[ a & 0x0000FFFF ] + bit_count16[ (a >> 32) & 0x0000FFFF ];
#endif

#ifndef TAIL_END_MASK
#define TAIL_END_MASK( x ) ~( 0xFFFFFFFFFFFFFFFC << ( 62 - ( x << 1 ) ) )
#endif

#elif PROCESSOR_WORD_SIZE == 32

#define U_HALF_MASK 0xAAAAAAAA
#define BLOCKS_PER_UNIT 2
#define TAIL_END 15

#ifndef AddToContingency
#define AddToContingency( n, t, a, b )                                      \
    t = (a & b);                                    \
    if( t ) {                                       \
        t = ((t | ((t & 0x55550000) >> 15)));   \
        n += bit_count16[ t & 0x0000FFFF ]; \
    }
#endif

#ifndef AddToDistribution
#define AddToDistribution(n, a )                                            \
    a = ((a | ((a & 0x55550000) >> 15)))            \
        n += bit_count16[ a & 0x0000FFFF ];
#endif

#ifndef TAIL_END_MASK
#define TAIL_END_MASK( x ) ( 0xFFFFFFFF >> ( 30 - ( x << 1 ) ) )
#endif

#elif PROCESSOR_WORD_SIZE == 16

#define U_HALF_MASK 0xAAAA
#define BLOCKS_PER_UNIT 1
#define TAIL_END 7

#ifndef AddToContingency
#define AddToContingency( n, t, a, b )                                      \
    t = (a & b);                                    \
    if( t ) {                                       \
        n += bit_count16[ t ];                      \
    }
#endif

#ifndef AddToDistribution
#define AddToDistribution(n, a )                                            \
    n += bit_count16[ a ];
#endif

#ifndef TAIL_END_MASK
#define TAIL_END_MASK( x ) ( 0xFFFC << ( 14 - ( x << 1 ) ) )
#endif

#else
#error "Invalid Processor Bit Width Defined"
#endif

static byte bit_count16[ 0x10000 ];
static joint_genotypes contingency_lookup[ 0x100000 ];
static byte skip_count[ 16 ];

/**
    The following Macro will split a bit string (VAL) of genotypes into the three separate bits strings (AA, AB, BB)
    where the even numbered bits are set IFF the genotype for corresponding genotype is set.

    Boolean logic below:

    For BB genotypes: take the high order bit AND it with the lower order bit. If they are they are both 1 (ie. bb= 11b),
    then lower order bit will be 1, and higher order bits will be cleared.

    Remember the XOR of the higher and lower order bits. If Genotype is  AA or AB, then the XOR will result in lower order bit being 1.
    If the Genotype is XX or BB, then the XOR will result in lower bit being 0.

    ANDing the XOR value with the input value will result in lower order bit being 1 IFF genotype is AA.

    Shifting the input value by 1 bit (dividing by 2) will result in the lower order bit being 1 IFF genotype is AB.
*/

inline void DecodeBitStrings( PWORD val, PWORD &aa, PWORD &ab, PWORD &bb) {
    PWORD half = ((val & U_HALF_MASK) >> 1);
    bb = ( val & half );
    PWORD tmp = (val ^ half);
    aa = val & tmp;
    ab = ((val >> 1) & tmp);
}

inline void IncrementFrequencyValue( frequency_table &ft, ulong aa, ulong ab, ulong bb) {
    aa = ((aa | ((aa & 0x5555000055550000) >> 15)));
    ft.aa += bit_count16[ aa & 0x0000FFFF ] + bit_count16[ (aa >> 32) & 0x0000FFFF ];
    ab = ((ab | ((ab & 0x5555000055550000) >> 15)));
    ft.ab += bit_count16[ ab & 0x0000FFFF ] + bit_count16[ (ab >> 32) & 0x0000FFFF ];
    bb = ((bb | ((bb & 0x5555000055550000) >> 15)));
    ft.bb += bit_count16[ bb & 0x0000FFFF ] + bit_count16[ (bb >> 32) & 0x0000FFFF ];
}

inline void IncrementFrequencyValue( frequency_table &ft, uint aa, uint ab, uint bb) {
    aa = ((aa | ((aa & 0x55550000) >> 15)));
    ft.aa += bit_count16[ aa & 0x0000FFFF ];
    ab = ((ab | ((ab & 0x55550000) >> 15)));
    ft.ab += bit_count16[ ab & 0x0000FFFF ];
    bb = ((bb | ((bb & 0x55550000) >> 15)));
    ft.bb += bit_count16[ bb & 0x0000FFFF ];
}

inline void IncrementFrequencyValue( frequency_table &ft, ushort aa, ushort ab, ushort bb) {
    ft.aa += bit_count16[ aa ];
    ft.ab += bit_count16[ ab ];
    ft.bb += bit_count16[ bb ];
}


class CompressedGenotypeTable3 : public GenoTable {
public:
    CompressedGenotypeTable3( indexer *markers, indexer *individs ) : GenoTable( markers, individs ), gt_lookup(NULL) {
        initialize();
    }

    DataBlock operator()( int r, int c );

    void addGenotype( int rIdx, int cIdx, const string &gt );
    void addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim );
    void addGenotypeRow( int rIdx, const char *p_begin, const char *p_end, char delim );

    ushort encodeGenotype( const string &gt );
    const char *decodeGenotype( ushort encoded_gt );

    bool isGenotypeHomozygous( ushort encoded_gt );

    void selectMarker( uint rIdx );

    void getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist );
    void getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet &ccs, CaseControlGenotypeDistribution &ccgd );

    void selectMarkerPair( uint maIdx, uint mbIdx );

    void getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct );
    void getContingencyTable( uint rIdx1, uint rIdx2, ushort *column_set, ContingencyTable &ct );
    void getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct );

    virtual ~CompressedGenotypeTable3();
protected:
    void initialize();
    void constructCountLookup();
    void constructContingencyLookup();


    char *gt_lookup, * err_lookup;
    uint gt_size, lookup_size;
    DataBlock **lookup;

    genotype_counts count_lookup[ 0x10000 ];
};

uint ones16( register uint x );


//    contingency_lookup = new joint_genotypes[ 0x100000 ]; // 16 * 256 * 256 == 2^4 * 2^8 * 2^8 == 2^20 == 0x100000
//    skip_count = new byte[ 16 ];

}
}

#endif // COMPRESSEDGENOTYPETABLE3_H
