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
#include "maf_func.h"

namespace libgwaspp {
namespace algorithms {

void computeMAF( GenoTable &gt ) {
    int _counts[ POSSIBLE_ENC ];

    int total = 0;
    int row_count = gt.row_size(), col_count = gt.column_size();
    for( int i = 0; i < row_count; ++i ) {
        memset( _counts, 0, POSSIBLE_ENC * sizeof( int ) );

        for( int j = 0; j < col_count; ++j ) {
            ++_counts[ gt( i, j )];
        }

        total = 0;
        for( int k = 0; k < POSSIBLE_ENC; ++k ) {
            if( _counts[ k ] > 0 ) {
                total += _counts[k];
            }
        }

        assert( total == col_count );
    }
}

void computeMAF( GeneticData *gd, const set<string> & marker_ids, const set<string> & individ_ids ) {
    set<int> r, c;

    int _counts[ POSSIBLE_ENC ];
    int total = 0, idx = 0;

    /**
    Convert Marker ids into table order row index
    */
    set<string>::const_iterator it = marker_ids.begin(), end = marker_ids.end();
    for( ; it != end; ++it ) {
        idx = gd->getGenotypedMarkerIndex( *it );
        assert( idx != -1 );

        r.insert( idx );
    }

    /**
    Covert Individual ids into table order column index
    */
    it = individ_ids.begin();
    end = individ_ids.end();
    for( ; it != end; ++it )  {
        idx = gd->getGenotypedIndividualIndex( *it );
        assert( idx != -1 );

        c.insert( idx );
    }

    /**
    All tables are only accessible via row and column indices
    */
    set<int>::iterator r_it = r.begin(), r_end = r.end();
    set<int>::iterator c_it, c_end;
    GenoTable &gt = *gd->getGenotypeTable();
    int individ_count = ( int ) individ_ids.size();
    for( ; r_it != r_end; ++r_it ) {
        memset( _counts, 0, POSSIBLE_ENC * sizeof( int ) );
        for( c_it = c.begin(), c_end = c.end() ; c_it != c_end; ++c_it ) {
            ++_counts[ gt(( *r_it ), ( *c_it ) )];
        }

        total = 0;
        for( int k = 0; k < POSSIBLE_ENC; ++k ) {
            if( _counts[ k ] > 0 ) {
                total += _counts[k];
            }
        }
        assert( total == individ_count );
    }
}

void maf( void *input, void *output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );

    int _counts[ POSSIBLE_ENC ];
    int total = 0;

    /**
    All tables are only accessible via row and column indices
    */
    set<int>::iterator r_it = inp.midx->begin(), r_end = inp.midx->end();
    set<int>::iterator c_it, c_end;
    GenoTable &gt = *inp.gd->getGenotypeTable();

    int individ_count = ( int ) inp.iidx->size();
    for( ; r_it != r_end; ++r_it ) {
        memset( _counts, 0, POSSIBLE_ENC * sizeof( int ) );
        for( c_it = inp.iidx->begin(), c_end = inp.iidx->end() ; c_it != c_end; ++c_it ) {
            ++_counts[ gt(( *r_it ), ( *c_it ) )];
        }

        total = 0;
        for( int k = 0; k < POSSIBLE_ENC; ++k ) {
            if( _counts[ k ] > 0 ) {
                total += _counts[k];
            }
        }
        assert( total == individ_count );
    }
}

void maf_all( void *input, void *output ) {

    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );

    int _counts[ POSSIBLE_ENC ];
    int total = 0, r_it, c_it, val;

    int marker_count = inp.gd->getGenotypedMarkersCount(), individ_count = inp.gd->getGenotypedIndividualsCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    ushort *p_begin = gt.possible_genotypes_begin(), * p_end = gt.possible_genotypes_end();

    cout << "List of possible Genotype Encodings: " << endl;
    for( ushort *c = p_begin, i = 0; c != p_end; c++, i++ ) {
        cout << i << "\t" << *c << " -> " << gt.decodeGenotype( *c )  << endl;
    }

    for( r_it = 0; r_it < marker_count; ++r_it ) {
        memset( _counts, 0, POSSIBLE_ENC * sizeof( int ) );
        for( c_it = 0; c_it < individ_count; ++c_it ) {
            ++_counts[ gt( r_it, c_it )];
        }

        total = 0;
        for( ushort *k = p_begin; k != p_end; k++ ) {
            if(( val = _counts[ *k ] ) > 0 ) {
                total += val;
            }
        }
        assert( total == individ_count );
    }
}

void maf_from_distribution( void *input, void *output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );

    int total = 0, r_it;

    int marker_count = inp.gd->getGenotypedMarkersCount(), individ_count = inp.gd->getGenotypedIndividualsCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    ushort *p_begin = gt.possible_genotypes_begin(), * p_end = gt.possible_genotypes_end();

    cout << "List of possible Genotype Encodings: " << endl;
    for( ushort *c = p_begin, i = 0; c != p_end; c++, i++ ) {
        cout << i << "\t" << *c << " -> " << gt.decodeGenotype( *c )  << endl;
    }

    frequency_table &dist = *gt.getGenotypeDistribution();
//    header_table &header = *gt.getGenotypeDistributionHeader();

//    INIT_LAPSE_TIME;

    for( r_it = 0; r_it < marker_count; ++r_it ) {
//        RECORD_START;
        gt.selectMarker( r_it );
//        RECORD_STOP;

        total = dist.xx;
        total += dist.aa;
        total += dist.ab;
        total += dist.bb;

//        cout << r_it <<  ". " << gt.decodeGenotype( header.xx ) << " - " << dist.xx << "; ";
//        cout << gt.decodeGenotype( header.aa ) << " - " << dist.aa << "; ";
//        cout << gt.decodeGenotype( header.ab ) << " - " << dist.ab << "; ";
//        cout << gt.decodeGenotype( header.bb ) << " - " << dist.bb << endl;
//
//        PRINT_LAPSE( "SELECT MARKER: " );

        assert( total == individ_count );
    }
}

void maf_all( IndexedInput &inp, void *output ) {
    int _counts[ POSSIBLE_ENC ];
    int total = 0, r_it, c_it;
    int val;
    int marker_count = inp.gd->getGenotypedMarkersCount(), individ_count = inp.gd->getGenotypedIndividualsCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    for( r_it = 0; r_it < marker_count; ++r_it ) {
        memset( _counts, 0, POSSIBLE_ENC * sizeof( int ) );
        for( c_it = 0; c_it < individ_count; ++c_it ) {
            ++_counts[ gt( r_it, c_it )];
        }

        total = 0;
        for( ushort *k = gt.possible_genotypes_begin(); k != gt.possible_genotypes_end(); k++ ) {
            if(( val = _counts[ *k ] ) > 0 ) {
                total += val;
            }
        }
        assert( total == individ_count );
    }
}

void compute_maf_perform( GeneticData *gd, ostream *out) {
    int marker_count = gd->getGenotypedMarkersCount();
    GenoTable &gt = *gd->getGenotypeTable();

    GenotypeDistribution dist;
    const frequency_table &ft = *dist.getDistribution();
    double tot, maf;

    for( int i = 0; i < marker_count; ++i ) {
        gt.getGenotypeDistribution( i, dist );

        tot = ft.aa;
        maf = 2.0 * tot;

        tot += ft.ab;
        maf += ft.ab;

        tot += ft.bb;

        maf /= tot;
        if( maf < 0.5 ) {
            maf = 1.0 - maf;
        }
    }
}

}
}
