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
#ifndef COMPRESSEDGENOTYPETABLE4_H
#define COMPRESSEDGENOTYPETABLE4_H

#include <iostream>
#include <cmath>

#include "common.h"
#include "util/index_set/indexer.h"
#include "genetics/genotype/geno_table.h"
#include "genetics/genotype/common_genotype.h"
#include "genetics/genotype/common_genotype_func.h"

using namespace std;
using namespace util;

namespace libgwaspp {
namespace genetics {

/**
    byte compressed genotype
    3 bits per Genotype
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

/**
 * Class: CompressedGenotypeTable4
 * Description: This class uses a 3-bit streaming approach to genotype compression
 * In other words, genotypes are divided into the 3 streams of bits. A stream
 * represents one of AA, AB, BB genotypes
 *
 * Testing claim that the streaming approach will allow for greater throughput
 */
class CompressedGenotypeTable4 : public GenoTable {
public:
    CompressedGenotypeTable4( indexer *markers, indexer *individs ) : GenoTable( markers, individs ), gt_lookup(NULL), m_cases(NULL), m_controls(NULL) {
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
    void getCaseControlGenotypeDistribution( uint rIdx, CaseControlGenotypeDistribution &ccgd );

    void selectMarkerPair( uint maIdx, uint mbIdx );
    void selectCaseControl( CaseControlSet &ccs );

    void getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct );
    void getContingencyTable( uint rIdx1, uint rIdx2, ushort *column_set, ContingencyTable &ct );
    void getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct );

    virtual ~CompressedGenotypeTable4();
protected:
    void initialize();
    void constructCountLookup();
    void constructContingencyLookup();


    char *gt_lookup, * err_lookup;
    uint gt_size, lookup_size;
    DataBlock **lookup;

    uint genotype_block_offset_ab, genotype_block_offset_bb;

    genotype_counts count_lookup[ 0x10000 ];

    DataBlock *m_cases, *m_controls;
    joint_genotypes contingency_lookup[ 0x100000 ];
    byte skip_count[ 16 ];
};

}
}

#endif // COMPRESSEDGENOTYPETABLE4_H
