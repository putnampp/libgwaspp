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

const char * g_Aheader[GENOTYPE_COUNT] = {"AA", "Aa", "aa", "xx" };
const char * g_Bheader[GENOTYPE_COUNT] = {"BB", "Bb", "bb", "xx" };

void printContingencyTable( const CONTIN_TABLE_T & ct, ostream & out, bool include_header ) {
    if( include_header ) {
        for( int j = 0; j < CONTIN_COLUMN_COUNT; ++j ) {
            out << "\t" << g_Bheader[ j ];
        }
        out << endl;
    }

    out << dec;
    for( int i = 0; i < CONTIN_ROW_COUNT; ++i ) {
        if( include_header )
            out << g_Aheader[ i ];
        for( int j = 0; j < CONTIN_COLUMN_COUNT; ++j ) {
            out << "\t" << ct.contin[ i * CONTIN_COLUMN_COUNT + j];
        }

        out << endl;
    }
}

void printFrequencyDistribution( const frequency_table & ft, ostream & out, bool include_header ) {
    if( include_header ) {
        for( int i = 0; i < GENOTYPE_COUNT; ++i ) {
            if( i )
                out << "\t";
            out << g_Aheader[i];
        }
        out << endl;
    }
    out << dec;
    for( int i = 0; i < GENOTYPE_COUNT; ++i ) {
        if( i )
            out << "\t";
        out << ft.freq[i];
    }
    out << endl;
}

ostream &operator<<(ostream &out, frequency_table & ft ) {
    out << dec;
    for( int i = 0; i < GENOTYPE_COUNT; ++i) {
        out << "\t" << ft.freq[i];
    }
    return out;
}

ostream &operator<<(ostream &out, CONTIN_TABLE_T & ct ) {
    /*out << dec << "\tBB\tBb\tbb" << endl;
    for( int i = 0; i < 9; ++i) {
        switch( i ) {
            case 0:
                out << "AA";
                break;
            case 3:
                out << "Aa";
                break;
            case 6:
                out << "aa";
                break;
            default:
                break;
        }

        out << "\t" << ct.contin[i];

        switch( i ) {
            case 2:
            case 5:
            case 8:
                out << endl;
            default:
                break;
        }
    }*/
    printContingencyTable( ct, out, true );
    return out;
}

ostream &operator<<(ostream &out, marginal_information &mi ) {
    //out << dec << "\tXX\tAA\tAB\tBB\t\tPbc_XX\tPbc_AA\tPbc_AB\tPbc_BB\t\tPca_XX\tPca_AA\tPca_AB\tPca_BB" << endl;
    //
    for( int i = 0; i < GENOTYPE_COUNT; ++i ) {
        out << "\t" << g_Aheader[i];
    }
    for( int i = 0; i < GENOTYPE_COUNT; ++i ) {
        out << "\tPbc_" << g_Aheader[i];
    }
    for( int i = 0; i < GENOTYPE_COUNT; ++i ) {
        out << "\tPca_" << g_Aheader[i];
    }
    out << endl;
    out << "Cases" << mi.cases;

    for( int i = 0; i < GENOTYPE_COUNT; ++i) {
        out << "\t" << mi.dPbc[i];
    }
    for( int i = 0; i < GENOTYPE_COUNT; ++i) {
        out << "\t" << mi.dPca[i];
    }
    out << endl;

    out << "Controls" << mi.controls;
    for( int i = GENOTYPE_COUNT; i < 2 * GENOTYPE_COUNT; ++i) {
        out << "\t" << mi.dPbc[i];
    }
    for( int i = GENOTYPE_COUNT; i < 2 * GENOTYPE_COUNT; ++i) {
        out << "\t" << mi.dPca[i];
    }
    out << endl;

    out << "Margins" << mi.margins << endl;

    out << "Entropy\t" << mi.dMarginalEntropy << endl;
    out << "Entropy_Y\t" << mi.dMarginalEntropy_Y << endl;
    return out;
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

int init_bit_count( byte * _count ) {
    byte * tmp = _count;
    for( uint i = 0; i < 0x10000; ++i, ++tmp ) {
        *tmp = ones16( i );
    }
    return 1;
}

void computeMarginalInformation( const frequency_table & _cases, const frequency_table & _ctrls, uint nIndivids, marginal_information & m) {
    CopyFrequencyTable( m.cases, _cases );
    CopyFrequencyTable( m.controls, _ctrls );

    uint nCases = m.cases.xx + m.cases.aa + m.cases.ab + m.cases.bb;
    uint nControls = m.controls.xx + m.controls.aa + m.controls.ab + m.controls.bb;

    assert( nCases + nControls == nIndivids );

    m.dMarginalEntropy = 0.0;
    m.dMarginalEntropy_Y = 0.0;

    uint * mar = &m.margins.freq[0];
    uint * _ca = &m.cases.freq[0];
    uint * _co = &m.controls.freq[0];

    double * pbc_ca = &m.dPbc[0];
    double * pbc_co = &m.dPbc[GENOTYPE_COUNT];
    double * pca_ca = &m.dPca[0];
    double * pca_co = &m.dPca[GENOTYPE_COUNT];

    double tmp;
    for( int i = 0; i < GENOTYPE_COUNT; ++i ) {
        *mar = *_ca + *_co;
        if( *mar > 0 ) {
            tmp = (double) *mar / (double) nIndivids;
            m.dMarginalEntropy += -(tmp) * log(tmp);
        }

        if( *_ca > 0 ) {
            tmp = (double) *_ca / nIndivids;
            m.dMarginalEntropy_Y += -(tmp) * log(tmp);
            *pbc_ca = (double) *_ca / (double) nCases;
            *pca_ca = (double) *_ca / *mar;
        }

        if( *_co > 0 ) {
            tmp = (double) *_co / nIndivids;
            m.dMarginalEntropy_Y += -(tmp) * log(tmp);

            *pbc_co = (double) *_co / (double) nControls;
            *pca_co = (double) *_co / (double) *mar;
        }

        ++mar; ++_ca; ++_co; ++pbc_ca; ++pbc_co; ++pca_ca; ++pca_co;
    }
}

}
}
