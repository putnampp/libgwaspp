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
#ifndef BASICGENOTYPETABLE_H
#define BASICGENOTYPETABLE_H

#include <cmath>
#include <cstring>

#include "common.h"
#include "genetics/genotype/geno_table.h"
#include "util/index_set/indexer.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

/**
    Basic Genotype Table
    No compression/compacted data
    1 - data block === 1 genotype

    /// Little-endian bit ordering
    /// BLOCK INDEX/GENOTYPE INDEX (ushort):          |       0       |       1       |  ...
    /// BASE INDEX (byte):                            |   0   |   1   |   2   |   3   |  ...
    /// BASE INDEX (bit):                             |15    8|7     0|15    8|7     0|  ...
*/

class BasicGenotypeTable : public GenoTable {
    public:
        BasicGenotypeTable(util::indexer * _rows, util::indexer * _columns) : GenoTable(_rows, _columns)  { initialize(); }

        DataBlock operator()( int r, int c);

        void addGenotype( int rIdx, int cIdx, const string &gt );
        void addGenotypeRow( int rIdx, string::const_iterator & it, string::const_iterator & it_end, char delim);
        void addGenotypeRow( int rIdx, const char * p_begin, const char * p_end, char delim );

        ushort encodeGenotype( const string &gt );
        const char *decodeGenotype( ushort encoded_gt );

        void selectMarker( uint rIdx );

        void getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist );
        void getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet & ccs, CaseControlGenotypeDistribution & ccgd );
        void getCaseControlGenotypeDistribution( uint rIdx, CaseControlGenotypeDistribution & ccgd ) { assert(false); }

        void selectMarkerPair( uint rIdx1, uint rIdx2 );

        void getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable & ct );
        void getContingencyTable( uint rIdx1, uint rIdx2, ushort *column_set, ContingencyTable &ct ) {
            assert( false );
        }
        void getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet & ccs, CaseControlContingencyTable & ccct );

        void selectCaseControl( CaseControlSet& ccs ) { assert( false ); }

        bool isGenotypeHomozygous( ushort enc );

        virtual ~BasicGenotypeTable();
    protected:
        void initialize();

        char * gt_lookup;
        ushort tmp_lookup[ 0x10000 ];
};

}
}

#endif // BASICGENOTYPETABLE_H
