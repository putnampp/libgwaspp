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
#ifndef COMPRESSEDGENOTYPETABLE_H
#define COMPRESSEDGENOTYPETABLE_H

#include <iostream>
#include <cmath>

#include "common.h"
#include "util/index_set/indexer.h"
#include "genetics/genotype/geno_table.h"
#include "genetics/genotype/common_genotype.h"


using namespace util;

namespace libgwaspp {
namespace genetics {

/**
    byte compressed genotype
    4-bits per base

    /// Little-endian bit ordering
    /// BLOCK INDEX (ushort):          |       0       |       1       |  ...
    /// DATA INDEX (byte):             |   0   |   1   |   2   |   3   |  ...
    /// GENOTYPE INDEX (half-byte):    | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 |  ...
    /// BASE INDEX (bit):              |15    8|7     0|15    8|7     0|  ...
*/

class CompressedGenotypeTable : public GenoTable {
    public:
        CompressedGenotypeTable( indexer *markers, indexer *individuals ) : GenoTable( markers, individuals ) { initialize();}

        DataBlock operator()( int r, int c );

        void addGenotype( int rIdx, int cIdx, const string &gt );
        void addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim );
        void addGenotypeRow( int rIdx, const char * p_begin, const char * p_end, char delim );

        ushort encodeGenotype( const string &gt );
        const char *decodeGenotype( ushort encoded_gt );

        bool isGenotypeHomozygous( ushort enc );

        void selectMarker( uint rIdx );

        void getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist );
        void getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet & ccs, CaseControlGenotypeDistribution & ccgd );

        void selectMarkerPair( uint rIdx1, uint rIdx2 );

        void getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable & ct );
        void getContingencyTable( uint rIdx1, uint rIdx2, ushort *column_set, ContingencyTable &ct ) {
            assert( false );
        }
        void getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet & ccs, CaseControlContingencyTable & ccct );

        virtual ~CompressedGenotypeTable();
    protected:
        void initialize();

        char * gt_lookup;

        union indices {
            ushort us;
            struct {
                byte hi, lo;
            };
        };
        indices tmp_lookup[ 0x10000 ];
    private:
};

}
}

#endif // COMPRESSEDGENOTYPETABLE_H
