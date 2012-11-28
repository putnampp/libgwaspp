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
#include "genetics/genotype/compressed_genotype_table2.h"

namespace libgwaspp {
namespace genetics {

void CompressedGenotypeTable2::initialize() {
    cout << "Initializing CompressedGenotypeTable2 ... " << endl;
    alphabet_size = ( int ) alphabet.length();
    possible_genotypes_size = pow( alphabet_size, ( double ) MAX_ALLELE_COUNT ) + 1;   // +1 for 0 -> "unknown" genotype

    cout << "Alphabet Size: " << alphabet_size << endl;
    cout << "Possible Genotype Size: " << possible_genotypes_size << endl;

    assert( possible_genotypes_size < ( 2 << BITS_PER_BLOCK ) );

    bits_per_data = 4;
    cout << "Bits per Data: " << bits_per_data << endl;

    assert(( BITS_PER_BLOCK % bits_per_data ) == 0 );   // assert that full data elements can be contained in a single block

    data_per_block = ( BITS_PER_BLOCK ) / bits_per_data;

    blocks_per_row = max_column / data_per_block + ((( max_column % data_per_block ) > 0 ) ? 1 : 0 ) + 1;
    total_block_count = blocks_per_row * max_row;
    bytes_per_row = blocks_per_row * BYTES_PER_BLOCK;

    data_size = total_block_count * BYTES_PER_BLOCK;

    cout << "Data per block: " << data_per_block << endl;
    cout << "Blocks per row: " << blocks_per_row << endl;
    cout << "Block count: " << total_block_count << endl;

    cout << "Total Table size: " << data_size << " (bytes)" << endl;

    if( data != NULL ) {
        delete [] data;
    }
    data = new DataBlock[ total_block_count ];
    memset( data, 0, data_size );

    if( gt_lookup != NULL ) {
        delete [] gt_lookup;
    }

    gt_lookup = new char[ MAX_ALLELE_COUNT + 1];
    gt_lookup[ MAX_ALLELE_COUNT ] = ( char )0;

    // assume every encoding is an error
    memset( transformations, ( byte )( alphabet_size ), 256 );

    int pattern = 0;
    for( string::const_iterator it = alphabet.begin(); it != alphabet.end(); ++it, pattern += 1 ) {
        transformations[( byte ) *it ] = pattern;
    }

    // build lookup table
    lookup = new DataBlock * [ alphabet_size + 1 ];
    beg = new ushort[ possible_genotypes_size];
    memset( beg, 0xFF, possible_genotypes_size * sizeof( ushort ) );
    end = beg + possible_genotypes_size;

    for( uint i = 0; i < alphabet_size + 1; ++i ) {
        lookup[ i ] = new DataBlock[ alphabet_size + 1 ];
        memset( lookup[i], 0xFF, ( alphabet_size + 1 ) * sizeof( DataBlock ) );
    }

    ushort enc = 0;
    for( uint i = 0, k = 0; i < alphabet_size; ++i ) {
        for( uint j = 0; j < alphabet_size; ++j ) {
            SetUshortAtDataBlock(lookup[i][j], enc);
            cout << GetUshortAtDataBlock(lookup[i][j]) << "\t";
            beg[k++] = enc++;
        }
        cout << endl;
    }

    gt_size =  MAX_ALLELE_COUNT + 1;
    int idx, offset;
    lookup_size = gt_size * possible_genotypes_size;

    gt_lookup = new char[ lookup_size ];    // allocate space for all possible genotypes as null-terminated character sequences
    err_lookup = gt_lookup + gt_size * ( possible_genotypes_size - 1 );
    memset( gt_lookup, ( char )0, lookup_size );
    for( uint i = 1, j = 0; i < possible_genotypes_size; ++i, ++j ) {
        idx = j;
        for( int k = 0, l = MAX_ALLELE_COUNT - 1; k < MAX_ALLELE_COUNT; ++k, --l ) {
            offset = idx % alphabet_size;
            gt_lookup[ j *gt_size + l ] = alphabet[offset];
            idx -= offset;
            idx /= alphabet_size;
        }
    }

    for( int k = 0; k < MAX_ALLELE_COUNT; ++k ) {
        err_lookup[ k ] = '0';
    }

    constructCountLookup();
    cout << "Initialized Genotype Lookups: " << lookup_size << " (bytes)" << endl;
}

DataBlock CompressedGenotypeTable2::operator()( int r, int c ) {
    ulong row_offset = ( ulong ) r * blocks_per_row;
    ulong block_idx = (( uint ) c >> 2 );

    ulong offset = row_offset + block_idx + 1;

    ushort val = GetUshortAtDataBlock(data[offset]);
    ushort header = GetUshortAtDataBlock(data[row_offset]);

    uint block_offset = ( c & 0x03 ) << 2;

    val &= ( 0xF000 >> block_offset );

    val >>= ( 12 - block_offset );
    DataBlock db;

    if( val < 0x0002 ) {
        if( val == 0x0001 ) {
            SetUshortAtDataBlock( db, (( header & 0x0F00 ) >> 8 ) );
        } else {
            SetUshortAtDataBlock( db, 0xFFFF );
        }
    } else {
        if( val == 0x0002 ) {
            SetUshortAtDataBlock( db, (( header & 0x00F0 ) >> 4 ) );
        } else {
            SetUshortAtDataBlock( db, ( header & 0x000F ) );
        }
    }
    return db;
}

void CompressedGenotypeTable2::addGenotype( int rIdx, int cIdx, const string &gt ) {
    ushort tmp_e = 0;
    ushort enc = encodeGenotype( gt );

    uint block_idx = ( cIdx >> 2 ) + 1; // Assume data_per_block == 8; + 1 => skip first block as it is the header block
    uint block_bit_offset = ( cIdx & 3 ) << 1;  // Assume data_per_block == 8

    uint row_offset = rIdx * blocks_per_row;

    DataBlock *tmp_data = data + row_offset;

    if( enc != 0xFFFF ) {
        DataBlock *tmp_head = tmp_data;
        ushort head_val = GetUshortAtDataBlockPtr(tmp_head);
        ushort geno_order = ( head_val & 0xF000 );
        ushort head_shift = 0, clear_code = 0;

        if( geno_order != 0x7000 ) {
            HeaderStateMachine( geno_order, enc, head_shift, clear_code, tmp_e, head_val );
            SetUshortAtDataBlockPtr( tmp_head, head_val );
        } else {
            if( isGenotypeHomozygous( enc ) ) {
                if(( head_val & 0x000F ) == enc ) {
                    tmp_e = 3;
                } else {
                    assert((( head_val & 0x0F00 ) >> 8 ) == enc );
                    tmp_e = 1;
                }
            } else {
                assert((( head_val & 0x00F0 ) >> 4 ) == enc );
                tmp_e = 2;
            }
        }
    } else {
        tmp_e = 0;
    }
    tmp_data += block_idx;

    GetUshortAtDataBlockPtr(tmp_data) &= ( 0xFFFF ^( 0x000F << block_bit_offset ) );     // clear current value at offset
    GetUshortAtDataBlockPtr(tmp_data) |= ( tmp_e << block_bit_offset );
}

void CompressedGenotypeTable2::addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim ) {
#if MAX_ALLELE_COUNT == 2
    ushort enc = 0, tmp_e = 0;
    ushort c1, c2;

    ushort enc_set[ 256 ];
    memset( enc_set, 0xFF, 256 * sizeof( ushort ) );

    DataBlock *tmp_data = data + rIdx * blocks_per_row;
    DataBlock *tmp_header = tmp_data;

    // make sure the data row is "unknown"
    memset( tmp_data, 0, bytes_per_row );

    // method offers fewer pointer dereferences
    ushort val = 0, shift = 12;

    ushort head_val = 0, geno_code = 0x0000;
    ushort clear_code = 0, head_shift = 0;

    ++tmp_data;
    do {
        if( shift > 12 ) {
            SetUshortAtDataBlockPtr( tmp_data, val);
            ++tmp_data;
            val = 0;
            shift = 12;
        }

        c1 = transformations[( byte )*it++];
        c2 = transformations[( byte )*it++];
        enc = GetUshortAtDataBlock(lookup[ c1 ][ c2 ]);

        // only modify columns which are "known"
        if( enc != 0xFFFF ) {
            if(( tmp_e = enc_set[ enc ] ) == 0xFFFF ) {
                assert( geno_code < 0x7000 );

                HeaderStateMachine( geno_code, enc, head_shift, tmp_e, clear_code, head_val );
                enc_set[ enc ] = tmp_e;
            }

            assert( tmp_e != 0 );
            val |= ( tmp_e << shift );
            tmp_e = 0;
        }

        shift -= 4;
    } while( it++ != it_end );

    SetUshortAtDataBlockPtr( tmp_data, val);

    // set header
    SetUshortAtDataBlockPtr( tmp_header, head_val);
#else
#error Incomplete implementation of adding genotype by row
#endif
}

void CompressedGenotypeTable2::addGenotypeRow( int rIdx, const char * p_begin, const char * p_end, char delim ) {
    assert(false);
}

void CompressedGenotypeTable2::selectMarker( uint rIdx ) {
    if( rIdx != current_dist_rIdx ) {
        resetSelectedMarker();

        current_dist_rIdx = rIdx;
        DataBlock *tmp_data = data + current_dist_rIdx * blocks_per_row;
        ushort header_val = GetUshortAtDataBlockPtr( tmp_data ), tail_end;
        ushort geno_code = ( header_val & 0xF000 );
        genotype_counts val;

        ParseHeader( gt_header, geno_code, header_val );

        ++tmp_data;
        for( uint i = 2; i < blocks_per_row; ++i, ++tmp_data ) {
            val.ui = count_lookup[ GetUshortAtDataBlockPtr(tmp_data) ].ui;
            gt_dist.xx += val.c0;
            gt_dist.aa += val.c1;
            gt_dist.ab += val.c2;
            gt_dist.bb += val.c3;
        }

        if(( tail_end = ( max_column & 3 ) ) == 3 ) {
            val.ui = count_lookup[ GetUshortAtDataBlockPtr(tmp_data) ].ui;
            gt_dist.xx += val.c0;
            gt_dist.aa += val.c1;
            gt_dist.ab += val.c2;
            gt_dist.bb += val.c3;
        } else {
            ushort v = GetUshortAtDataBlockPtr(tmp_data);
            for( uint i = 0; i < tail_end; ++i, v <<= 4 ) {
                ++gt_dist.freq[( v & 0xF000 ) >> 12 ];
            }
        }
        row_selected = true;
    }
}

void CompressedGenotypeTable2::selectMarkerPair( uint rIdx1, uint rIdx2 ) {
    if( !(( maIdx == rIdx1 && mbIdx == rIdx2 ) || ( maIdx == rIdx2 && mbIdx == rIdx1 ) ) ) {
        resetSelectedMarkerPair();

        maIdx = rIdx1;
        mbIdx = rIdx2;

        DataBlock *ma_tmp_data = data + maIdx * blocks_per_row;
        DataBlock *mb_tmp_data = data + mbIdx * blocks_per_row;

        ushort ma_header_val = GetUshortAtDataBlockPtr( ma_tmp_data ), mb_header_val = GetUshortAtDataBlockPtr( mb_tmp_data );
        int tail_end = 0;
        ushort ma_geno_code = ( ma_header_val & 0xF000 ), mb_geno_code = ( mb_header_val & 0xF000 );
        joint_genotypes val;

        ParseHeader( ma_header, ma_geno_code, ma_header_val );
        ParseHeader( mb_header, mb_geno_code, mb_header_val );

        ++ma_tmp_data; ++mb_tmp_data;   // move past header
        uint base_offset = 0x30000;
        ushort a_val, b_val;
        for( uint i = 2; i < blocks_per_row; ++i, ++ma_tmp_data, ++mb_tmp_data ) {
            a_val = *ma_tmp_data;
            b_val = *mb_tmp_data;

            val = contingency_lookup[ base_offset | ( a_val & 0xFF00 ) | (( b_val & 0xFF00 ) >> 8 )];
            ContingencyAddJointGenotype( _contingency, val );

            val = contingency_lookup[ base_offset | (( a_val & 0x00FF ) << 8 ) | ( b_val & 0x00FF )];
            ContingencyAddJointGenotype( _contingency, val );
        }

        tail_end = (max_column & 3);

        a_val = *ma_tmp_data;
        b_val = *mb_tmp_data;

        if( tail_end >= 2 ) {
            val = contingency_lookup[ base_offset | ( a_val & 0xFF00 ) | (( b_val & 0xFF00 ) >> 8 )];

            ContingencyAddJointGenotype( _contingency, val );
            tail_end -= 2;
        } else {
            a_val >>= 8;
            b_val >>= 8;
        }


        if ( tail_end > 0 ) {
            if( tail_end == 1 ) {
                base_offset = 0x20000;
            } else {
                assert( tail_end == 2 );
                base_offset = 0x30000;
            }

            tail_end = 2 - tail_end;
            val = contingency_lookup[ base_offset | (( a_val & 0x00FF ) << 8 ) | ( b_val & 0x00FF )];

            ContingencyAddJointGenotype( _contingency, val );
            _contingency.n9 -= tail_end;
        }
    }
}

void CompressedGenotypeTable2::getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist ) {
    assert( false );
}

void CompressedGenotypeTable2::getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet &ccs, CaseControlGenotypeDistribution &ccgd ) {
    assert( false );
}

void CompressedGenotypeTable2::getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct ) {
    assert( false );
}
void CompressedGenotypeTable2::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct ) {
    assert( false );
}
void CompressedGenotypeTable2::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlContingencyTable &ccct ) {
    assert(false);
}

void CompressedGenotypeTable2::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, const marginal_information &m1, const marginal_information &m2, CaseControlContingencyTable & ccct ) {
    assert(false);
}

bool CompressedGenotypeTable2::isGenotypeHomozygous( ushort enc ) {
    // AA == 0; CC == 5; GG == 10; TT == 15
    return (( enc != 0xFFFF ) && (( enc == 0 ) || ( enc == 5 ) || ( enc == 10 ) || ( enc == 15 ) ) );
}

ushort CompressedGenotypeTable2::encodeGenotype( const string &gt ) {
    assert( gt.length() == MAX_ALLELE_COUNT );

    return (ushort) lookup[ transformations[( byte ) gt[0] ] ][ transformations[( byte ) gt[1] ] ];
}

const char *CompressedGenotypeTable2::decodeGenotype( ushort encoded_gt ) {
    if( encoded_gt == 0xFFFF ) {
        return err_lookup;
    }
    return gt_lookup + encoded_gt * gt_size;
}

void CompressedGenotypeTable2::constructCountLookup( ) {
    genotype_counts counts;
    for( uint i = 0; i < 0x10000; ++i ) {
        counts.ui = 0;

        for( uint j = 0, k = i; j < 4; ++j, k >>= 4 ) {
            switch( k & 0x000F ) {
                case 0x0003:
                    ++counts.c3;
                    break;
                case 0x0002:
                    ++counts.c2;
                    break;
                case 0x0001:
                    ++counts.c1;
                    break;
                case 0x0000:
                    ++counts.c0;
                    break;
                default:
                    break;
            }
        }

        count_lookup[i].ui = counts.ui;
    }
}

void CompressedGenotypeTable2::constructContingencyLookup() {
    joint_genotypes jg;

    memset( contingency_lookup, 0, 4 * 256 * 256 * sizeof( joint_genotypes ) );
    byte ma, mb;

    uint base_offset, offset;
    for( uint k = 1, mask = 0x00; k < 4; ++k ) {
        mask = 0;
        mask |= (( k & 0x01 ) ? 0xF0 : 0 );
        mask |= (( k & 0x02 ) ? 0x0F : 0 );

        for( uint i = 0; i < 256; ++i ) {
            base_offset = (( k << 16 ) | ( i << 8 ) );
            for( uint j = 0; j < 256; ++j ) {
                // clear data
                ResetJointGenotype( jg );

                for( int bit_mask = 0xF0, shift = 4, idx = 0; shift >= 0; bit_mask >>= 4, shift -= 4, ++idx ) {

                    ma = (( i & ( bit_mask & mask ) ) >> shift );
                    mb = (( j & ( bit_mask & mask ) ) >> shift );

                    SetJointGenotype( jg, ma, mb, idx );
                }

                offset = base_offset | j;
                CopyJointGenotype( contingency_lookup[ offset ], jg );
            }
        }
    }
}

CompressedGenotypeTable2::~CompressedGenotypeTable2() {
    //dtor
    delete [] gt_lookup;

    for( uint i = 0; i < alphabet_size + 1; i++ ) {
        delete [] lookup[ i ];
    }

    delete [] lookup;
}

}
}



