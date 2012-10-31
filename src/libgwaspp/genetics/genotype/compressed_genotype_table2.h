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
#ifndef COMPRESSEDGENOTYPETABLE2_H
#define COMPRESSEDGENOTYPETABLE2_H

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
    4 bits per Genotype
    1 block per Record to identify which genotype code is used as "UNKNOWN" for record

    Little-endian bit ordering
    BLOCK INDEX (ushort):          | HEADER |       0       |       1       |  ...
    GENOTYPE INDEX (half-byte):    |        | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |  ...
    CODE INDEX (bit):              |15     0|15    8|7     0|15    8|7     0|  ...

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

    Unknown Genotype value = 0xFFFF
*/

class CompressedGenotypeTable2 : public GenoTable {
    public:
        CompressedGenotypeTable2( indexer *markers, indexer *individuals ) : GenoTable( markers, individuals ) { initialize();}

        DataBlock operator()( int r, int c );

        void addGenotype( int rIdx, int cIdx, const string &gt );
        void addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim );
        void addGenotypeRow( int rIdx, const char * p_begin, const char * p_end, char delim );

        ushort encodeGenotype( const string &gt );
        const char *decodeGenotype( ushort encoded_gt );

        bool isGenotypeHomozygous( ushort enc );

        void selectMarker( uint rIdx );

        void getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist );
        void getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet &ccs, CaseControlGenotypeDistribution &ccgd );

        void selectMarkerPair( uint rIdx1, uint rIdx2 );

        void getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct );
        void getContingencyTable( uint rIdx1, uint rIdx2, ushort *column_set, ContingencyTable &ct ) {
            assert( false );
        }
        void getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct );

        virtual ~CompressedGenotypeTable2();
    protected:
        void initialize();
        void constructCountLookup();
        void constructContingencyLookup();

//        void refactorUnknownOfRow( uint row );

        char *gt_lookup, * err_lookup;
        DataBlock **lookup;

        genotype_counts count_lookup[ 0x10000 ];
        joint_genotypes contingency_lookup[ 0x40000 ]; // 4 * 256 * 256 == 2^2 * 2^8 * 2^8 == 2^18 == 0x40000

        uint gt_size, lookup_size;
};

}
}

#endif // COMPRESSEDGENOTYPETABLE2_H
