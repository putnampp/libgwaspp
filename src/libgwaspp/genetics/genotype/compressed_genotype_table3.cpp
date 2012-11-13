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
#include "genetics/genotype/compressed_genotype_table3.h"

namespace libgwaspp {
namespace genetics {

void CompressedGenotypeTable3::initialize() {
    cout << "Initializing CompressedGenotypeTable3 ... " << endl;
    alphabet_size = ( int ) alphabet.length();
    possible_genotypes_size = pow( alphabet_size, ( double ) MAX_ALLELE_COUNT ) + 1;   // +1 -> "unknown" genotype

    cout << "Alphabet Size: " << alphabet_size << endl;
    cout << "Possible Genotype Size: " << possible_genotypes_size << endl;

    assert( (int)possible_genotypes_size < ( 2 << BITS_PER_BLOCK ) );

    bits_per_data = 2;
    cout << "Bits per Data: " << bits_per_data << endl;

    assert(( BITS_PER_BLOCK % bits_per_data ) == 0 );   // assert that full data elements can be contained in a single block

    data_per_block = ( BITS_PER_BLOCK ) / bits_per_data;

    // always pad the blocks per row by 1.
    // safe assumption that max_columns will not likely be evenly divisible by data_per_block 
    blocks_per_row = max_column / data_per_block + 1;
    int block_bit_width = (sizeof(DataBlock) << 3);
    assert( (PROCESSOR_WORD_SIZE % block_bit_width) == 0 );

    // further pad the number of data blocks per row for processor word size efficiency
    int block_per_pword = (PROCESSOR_WORD_SIZE / block_bit_width);
    if( block_per_pword > 1 && blocks_per_row % block_per_pword ) {
        blocks_per_row += ( block_per_pword - (blocks_per_row % block_per_pword));
    }

    // finally, add an additional block for genotype headers
    ++blocks_per_row;
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
    data_size = max_row * bytes_per_row;

    cout << "Total Table size: " << data_size << " (bytes)" << endl;

    data = new DataBlock[ data_size ];
    memset( data, 0, data_size );

    if( gt_lookup != NULL ) {
        delete [] gt_lookup;
    }

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
            SetUshortAtDataBlock( lookup[i][j], enc );
            beg[k++] = enc++;
        }
    }

    for( uint i = 0; i < alphabet_size + 1; ++i ) {
        for( uint j = 0; j < alphabet_size + 1; ++j ) {
            cout << GetUshortAtDataBlock( lookup[i][j] ) << "\t";
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
    constructContingencyLookup();

    cout << "Initialized Genotype Lookups: " << lookup_size << " (bytes)" << endl;
}

void CompressedGenotypeTable3::addGenotype( int rIdx, int cIdx, const string &gt ) {
    ushort tmp_e = 0;
    ushort enc = encodeGenotype( gt );

    uint block_idx = ( cIdx >> 3 ) + 1; // Assume data_per_block == 8; + 1 => skip first block as it is the header block
    uint block_bit_offset = (( cIdx & 7 ) << 1);  // Assume data_per_block == 8

    uint row_offset = rIdx * blocks_per_row;

    DataBlock *tmp_data = data + row_offset;

    if( enc != 0xFFFF ) {
        DataBlock *tmp_head = tmp_data;
        ushort head_val = GetUshortAtDataBlockPtr( tmp_head );
        ushort geno_order = ( head_val & 0xF000 );
        ushort head_shift = 0, clear_code = 0;

        if( geno_order != 0x7000 ) {
            HeaderStateMachine( geno_order, enc, head_shift, clear_code, tmp_e, head_val );
//            *tmp_head = head_val;   // set header
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

    GetUshortAtDataBlockPtr( tmp_data ) &= ( 0xFFFF ^ ( 0x0003 << block_bit_offset ) );  // clear current value at offset
    GetUshortAtDataBlockPtr( tmp_data ) |= ( tmp_e << block_bit_offset );
}

void CompressedGenotypeTable3::addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim ) {
#if MAX_ALLELE_COUNT == 2
    if( it >= it_end ) return;

    ushort enc = 0, tmp_e = 0;
    ushort c1, c2;

    ushort enc_set[ 256 ];
    memset( enc_set, 0xFF, 256 * sizeof( ushort ) );

    DataBlock *tmp_data = data + rIdx * blocks_per_row;
    DataBlock *tmp_header = tmp_data;

    // make sure the data row is "unknown"
    memset( tmp_data, 0, bytes_per_row );

    // method offers fewer pointer dereferences
    ushort val = 0, shift = 0;

    ushort head_val = 0, geno_code = 0x0000;
    ushort clear_code = 0, head_shift = 0;

    ++tmp_data;         // skip header
    do {
        if( shift > 14 ) {
//            *tmp_data = val;
            SetUshortAtDataBlockPtr( tmp_data, val );
            ++tmp_data;
            val = 0;
            shift = 0;
        }

        c1 = transformations[( byte )*it++];
        c2 = transformations[( byte )*it++];

        enc = GetUshortAtDataBlock( lookup[ c1 ][ c2 ] );

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

        shift += 2;
    } while( ++it < it_end );

//    *tmp_data = val;
    SetUshortAtDataBlockPtr( tmp_data, val );

    // set header
//    *tmp_header = head_val;
    SetUshortAtDataBlockPtr( tmp_header, head_val );
#else
#error "Incomplete implementation of adding genotype by row"
#endif
}

void CompressedGenotypeTable3::addGenotypeRow( int rIdx, const char *p_begin, const char *p_end, char delim ) {
#if MAX_ALLELE_COUNT == 2
    if( p_begin >= p_end) return;

    ushort enc = 0, tmp_e = 0;
    ushort c1, c2;

    ushort enc_set[ 256 ];
    memset( enc_set, 0xFF, 256 * sizeof( ushort ) );

    DataBlock *tmp_data = data + rIdx * blocks_per_row;
    DataBlock *tmp_header = tmp_data;

    // make sure the data row is "unknown"
    memset( tmp_data, 0, bytes_per_row );

    // method offers fewer pointer dereferences
    ushort val = 0, shift = 0;

    ushort head_val = 0, geno_code = 0x0000;
    ushort clear_code = 0, head_shift = 0;

    ++tmp_data;
    const char *tmp_p = p_begin;
    do {
        if( shift > 14 ) {
            SetUshortAtDataBlockPtr( tmp_data, val );
            ++tmp_data;
            val = 0;
            shift = 0;
        }

        c1 = transformations[( byte )*tmp_p++];
        c2 = transformations[( byte )*tmp_p++];

        enc = GetUshortAtDataBlock( lookup[ c1 ][ c2 ] );

        // only modify columns which are "known"
        if( enc != 0xFFFF ) {
            if(( tmp_e = enc_set[ enc ] ) == 0xFFFF ) {
                assert( geno_code < 0x7000 );

                HeaderStateMachine( geno_code, enc, head_shift, tmp_e, clear_code, head_val );

                enc_set[ enc ] = tmp_e;
            }

            assert( tmp_e != 0 );
            //cout << "(" << tmp_e << ")";
            val |= ( tmp_e << shift );
            tmp_e = 0;
        }
        //cout << endl;
        shift += 2;
    } while( ++tmp_p < p_end );

    SetUshortAtDataBlockPtr( tmp_data, val );
    // set header
    SetUshortAtDataBlockPtr( tmp_header, head_val );

    //cout << hex << "Header: " <<  head_val << endl;
    //cout << dec;
#else
#error "Incomplete implementation of adding genotype by row"
#endif
}

ushort CompressedGenotypeTable3::encodeGenotype( const string &gt ) {
#if MAX_ALLELE_COUNT == 2
    assert( gt.length() == MAX_ALLELE_COUNT );

    return ( ushort ) lookup[ transformations[( byte ) gt[0] ] ][ transformations[( byte ) gt[1] ] ];
#else
    int idx_forms[ MAX_ALLELE_COUNT ];
    int found_count = 0, idx = 0;
    char c;
    for( string::iterator it = alphabet.begin(); it != alphabet.end() && found_count < MAX_ALLELE_COUNT; ++it, ++idx ) {
        c = *it;
        for( int i = 0; i < MAX_ALLELE_COUNT && found_count < MAX_ALLELE_COUNT; ++i ) {
            if( gt[i] == c ) {
                idx_forms[i] = idx;
                ++found_count;
            }
        }
    }

    if( found_count != MAX_ALLELE_COUNT ) return 0;

    idx = 0;
    for( int i = 0; i < MAX_ALLELE_COUNT && idx <= 255; ++i ) {
        idx *= alphabet_size;
        idx += idx_forms[i];
    }

    ++idx;
    if( idx > 255 ) return 0;
    return idx;
#endif
}

DataBlock CompressedGenotypeTable3::operator()( int r, int c ) {
    ulong row_offset = ( ulong ) r * blocks_per_row;
    ulong block_idx = (( uint ) c >> 3 );

    ulong offset = row_offset + block_idx + 1;

    ushort val = GetUshortAtDataBlock( data[offset] );
    ushort header = GetUshortAtDataBlock( data[row_offset] );

    uint block_offset = ( c & 0x07 ) << 1;

    //val &= ( 0x0003 << block_offset );
    val = ((val >> block_offset) & 0x0003);

    DataBlock db;

    /*if( val < 0x0002 ) {
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
    }*/
    switch( val ) {
    case 3:
        SetUshortAtDataBlock( db, (header & 0x000F) );
        break;
    case 2:
        SetUshortAtDataBlock( db, ((header & 0x00F0) >> 4) );
        break;
    case 1:
        SetUshortAtDataBlock( db, ((header & 0x0F00) >> 8 ) );
        break;
    default:
        SetUshortAtDataBlock( db, 0xFFFF );
        break;
    }
    return db;
}


void CompressedGenotypeTable3::selectMarker( uint rIdx ) {
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
            val = count_lookup[ GetUshortAtDataBlockPtr( tmp_data )];
            gt_dist.xx += val.c0;
            gt_dist.aa += val.c1;
            gt_dist.ab += val.c2;
            gt_dist.bb += val.c3;
        }

        if(( tail_end = ( max_column & 7 ) ) == 7 ) {
            val = count_lookup[ GetUshortAtDataBlockPtr( tmp_data )];
            gt_dist.xx += val.c0;
            gt_dist.aa += val.c1;
            gt_dist.ab += val.c2;
            gt_dist.bb += val.c3;
        } else {
            ushort v = GetUshortAtDataBlockPtr( tmp_data );
            for( uint i = 0; i < tail_end; ++i, v <<= 2 ) {
                ++gt_dist.freq[( v & 0xC000 ) >> 14 ];
            }
        }
        row_selected = true;
    }
}


void CompressedGenotypeTable3::selectMarkerPair( uint rIdx1, uint rIdx2 ) {
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

        ++ma_tmp_data;
        ++mb_tmp_data;   // move past header
        uint base_offset = 0xF0000;
        ushort a_val, b_val;
        for( uint i = 2; i < blocks_per_row; ++i, ++ma_tmp_data, ++mb_tmp_data ) {
            a_val = *ma_tmp_data;
            b_val = *mb_tmp_data;

            val = contingency_lookup[ base_offset | ( a_val & 0xFF00 ) | (( b_val & 0xFF00 ) >> 8 )];
            ContingencyAddJointGenotype( _contingency, val );

            val = contingency_lookup[ base_offset | (( a_val & 0x00FF ) << 8 ) | ( b_val & 0x00FF )];
            ContingencyAddJointGenotype( _contingency, val );
        }

        // tail_end is the first index which is empty
        tail_end = ( max_column & 7 );

        a_val = *ma_tmp_data;
        b_val = *mb_tmp_data;

        if( tail_end >= 4 ) {
            val = contingency_lookup[ base_offset | ( a_val & 0xFF00 ) | (( b_val & 0xFF00 ) >> 8 )];

            ContingencyAddJointGenotype( _contingency, val );

            tail_end -= 4;
        } else {
            a_val >>= 8;
            b_val >>= 8;
        }


        if( tail_end > 0 ) {
            if( tail_end == 1 ) {
                base_offset = 0x80000;
            } else if( tail_end == 2 ) {
                base_offset = 0xC0000;
            } else if( tail_end == 3 ) {
                base_offset = 0xE0000;
            } else {
                assert( tail_end == 4 );
                base_offset = 0xF0000;
            }

            tail_end = 4 - tail_end;    // padding genotypes

            val = contingency_lookup[ base_offset | (( a_val & 0x00FF ) << 8 ) | ( b_val & 0x00FF )];

            ContingencyAddJointGenotype( _contingency, val );
            _contingency.n9 -= tail_end;
        }
    }
}


void CompressedGenotypeTable3::selectCaseControl( CaseControlSet & ccs ) {
    if ( m_cases == NULL ) {
        m_cases = new DataBlock[data_size]; 
    }

    if ( m_controls == NULL ) {
        m_controls = new DataBlock[data_size]; 
    }

    // Copy contents of data into case and controls
    memcpy( m_cases, data, data_size );
    memcpy( m_controls, data, data_size );

    // Mask out cases and controls
    const PWORD *case_ptr = reinterpret_cast< const PWORD * >( ccs.case_begin() );
    const PWORD *ctrl_ptr = reinterpret_cast< const PWORD * >( ccs.control_begin() );

    const PWORD *tmp_case, *tmp_ctrl;

    // for every row
    for( uint i = 0, offset = 1; i < (uint)max_row; ++i, offset += blocks_per_row ) {
        // locate the start of the data segment for each row 
        PWORD *case_data = reinterpret_cast< PWORD * >( m_cases + offset );
        PWORD *ctrl_data = reinterpret_cast< PWORD * >( m_controls + offset );

        // reset the mask pointers
        tmp_case = case_ptr;
        tmp_ctrl = ctrl_ptr;
        // for every data block
        for( uint j = 1; j < blocks_per_row; j += BLOCKS_PER_PWORD ) {
            // mask out all unnecessary data columns
            *case_data++ &= *tmp_case++;
            *ctrl_data++ &= *tmp_ctrl++;
        } 
    }
}

void CompressedGenotypeTable3::getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist ) {
    PWORD *tmp_data = reinterpret_cast< PWORD * >(data + rIdx * blocks_per_row + 1);
    PWORD val, _aa, _ab, _bb, uhalf, lhalf;

    frequency_table joint_gt;
    ResetFrequencyTable(joint_gt);

    for( uint i = 1; i < blocks_per_row; i += BLOCKS_PER_PWORD ) {
        val = *tmp_data++;

        /*cout << hex;
        cout.width(16);
        cout << val;*/
        DecodeBitStrings2BitBlock(val, _aa, _ab, _bb);
        /*cout << hex << "\t";
        cout.width(16);
        cout << _aa << "\t";
        cout.width(16);
        cout << _ab << "\t";
        cout.width(16);
        cout << _bb;*/
        IncrementFrequencyValue2BitBlock(joint_gt, _aa, _ab, _bb);
        //cout << dec << "\t" << joint_gt.aa << "\t" << joint_gt.ab << "\t" << joint_gt.bb << endl;
    }

    /*
     * added processor width padding to data blocks per row
     * as a result ther should no longer be any tails
    PWORD tail_end = ( max_column & TAIL_END );
    if( tail_end > 0 ) {
        val = *tmp_data++;

        DecodeBitStrings2BitBlock(val, _aa, _ab, _bb);
        IncrementFrequencyValue(joint_gt, _aa, _ab, _bb);
    }
    */

    dist.setDistribution( joint_gt );
}

void CompressedGenotypeTable3::getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet &ccs, CaseControlGenotypeDistribution &ccgd ) {
    PWORD *tmp_data = reinterpret_cast< PWORD * >(data + rIdx * blocks_per_row + 1);

    PWORD val, case_val, ctrl_val;
    register PWORD _aa, _bb, _ab;

    frequency_table case_gt, ctrl_gt;
    ResetFrequencyTable( case_gt );
    ResetFrequencyTable( ctrl_gt );

    const PWORD *case_ptr = reinterpret_cast< const PWORD * >( ccs.case_begin() );
    const PWORD *ctrl_ptr = reinterpret_cast< const PWORD * >( ccs.control_begin() );

    // slip header block
    for( uint i = 1; i < blocks_per_row; i += BLOCKS_PER_PWORD ) {
        val = *tmp_data++;
        case_val = val & *case_ptr++;
        ctrl_val = val & *ctrl_ptr++;

        DecodeBitStrings2BitBlock( case_val, _aa, _ab, _bb );
        IncrementFrequencyValue2BitBlock(case_gt, _aa, _ab, _bb);

        DecodeBitStrings2BitBlock(ctrl_val, _aa, _ab, _bb );
        IncrementFrequencyValue2BitBlock(ctrl_gt, _aa, _ab, _bb);
    }

/*
    PWORD tail_end = ( max_column & TAIL_END );

    if( tail_end > 0 ) {
        val = *tmp_data++;

        case_val = val & *case_ptr++;
        ctrl_val = val & *ctrl_ptr++;

        DecodeBitStrings2BitBlock( case_val, _aa, _ab, _bb );
        IncrementFrequencyValue(case_gt, _aa, _ab, _bb);

        DecodeBitStrings2BitBlock( ctrl_val, _aa, _ab, _bb );
        IncrementFrequencyValue(ctrl_gt, _aa, _ab, _bb);
    }
*/
    ccgd.setCaseDistribution(case_gt);
    ccgd.setControlDistribution(ctrl_gt);
}

void CompressedGenotypeTable3::getCaseControlGenotypeDistribution( uint rIdx, CaseControlGenotypeDistribution &ccgd ) {
    ulong offset = rIdx * blocks_per_row + 1;
    PWORD *tmp_case_data = reinterpret_cast< PWORD * >( m_cases + offset );
    PWORD *tmp_ctrl_data = reinterpret_cast< PWORD * >( m_controls + offset );

    PWORD case_val, ctrl_val;
    PWORD _aa, _ab, _bb, uhalf, lhalf;

    frequency_table case_gt, ctrl_gt;
    ResetFrequencyTable( case_gt );
    ResetFrequencyTable( ctrl_gt );

    for( uint i = 1; i < blocks_per_row; i += BLOCKS_PER_PWORD ) {
        case_val = *tmp_case_data++;
        ctrl_val = *tmp_ctrl_data++;

        DecodeBitStrings2BitBlock(case_val, _aa, _ab, _bb );
        IncrementFrequencyValue2BitBlock(case_gt, _aa, _ab, _bb );

        DecodeBitStrings2BitBlock(ctrl_val, _aa, _ab, _bb );
        IncrementFrequencyValue2BitBlock( ctrl_gt, _aa, _ab, _bb);
    }

    ccgd.setCaseDistribution(case_gt);
    ccgd.setControlDistribution(ctrl_gt);
}


void CompressedGenotypeTable3::getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct ) {
    PWORD *ma_tmp_data = reinterpret_cast< PWORD * >( data + rIdx1 * blocks_per_row + BLOCKS_PER_PWORD );
    PWORD *mb_tmp_data = reinterpret_cast< PWORD * >( data + rIdx2 * blocks_per_row + BLOCKS_PER_PWORD );

    ct.setMarkerAIndex( rIdx1 );
    ct.setMarkerBIndex( rIdx2 );

    contingency_table contingency;
    ResetContingencyTable( contingency );

    register PWORD a_val, b_val;
    register PWORD a_aa, a_bb, b_aa, b_bb, a_ab, b_ab;
    PWORD tmp_val, uhalf, lhalf;

    for( uint i = 1; i < blocks_per_row; i += BLOCKS_PER_PWORD ) {
        a_val = *ma_tmp_data++;
        b_val = *mb_tmp_data++;

        DecodeBitStrings2BitBlock( a_val, a_aa, a_ab, a_bb );
        DecodeBitStrings2BitBlock( b_val, b_aa, b_ab, b_bb );

        AddToContingency2BitBlock( contingency.n0, tmp_val, a_aa, b_aa );
        AddToContingency2BitBlock( contingency.n1, tmp_val, a_aa, b_ab );
        AddToContingency2BitBlock( contingency.n2, tmp_val, a_aa, b_bb );
        AddToContingency2BitBlock( contingency.n3, tmp_val, a_ab, b_aa );
        AddToContingency2BitBlock( contingency.n4, tmp_val, a_ab, b_ab );
        AddToContingency2BitBlock( contingency.n5, tmp_val, a_ab, b_bb );
        AddToContingency2BitBlock( contingency.n6, tmp_val, a_bb, b_aa );
        AddToContingency2BitBlock( contingency.n7, tmp_val, a_bb, b_ab );
        AddToContingency2BitBlock( contingency.n8, tmp_val, a_bb, b_bb );
    }
/*
    PWORD tail_end = ( max_column & TAIL_END );
    if( tail_end > 0 ) {
        tmp_val = ( 0xFFFFFFFFFFFFFFFC << ( 62 - ( tail_end << 1 ) ) );

        a_val = *ma_tmp_data++;
        a_val &= tmp_val;

        b_val = *mb_tmp_data++;
        b_val &= tmp_val;

        DecodeBitStrings2BitBlock( a_val, a_aa, a_ab, a_bb );
        DecodeBitStrings2BitBlock( b_val, b_aa, b_ab, b_bb );

        AddToContingency2BitBlock( contingency.n0, tmp_val, a_aa, b_aa );
        AddToContingency2BitBlock( contingency.n1, tmp_val, a_aa, b_ab );
        AddToContingency2BitBlock( contingency.n2, tmp_val, a_aa, b_bb );
        AddToContingency2BitBlock( contingency.n3, tmp_val, a_ab, b_aa );
        AddToContingency2BitBlock( contingency.n4, tmp_val, a_ab, b_ab );
        AddToContingency2BitBlock( contingency.n5, tmp_val, a_ab, b_bb );
        AddToContingency2BitBlock( contingency.n6, tmp_val, a_bb, b_aa );
        AddToContingency2BitBlock( contingency.n7, tmp_val, a_bb, b_ab );
        AddToContingency2BitBlock( contingency.n8, tmp_val, a_bb, b_bb );
    }
*/
    ct.setContingency( contingency );
}

void CompressedGenotypeTable3::getContingencyTable( uint rIdx1, uint rIdx2, ushort *column_set, ContingencyTable &ct ) {

}

void CompressedGenotypeTable3::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct ) {

    PWORD *ma_tmp_data = reinterpret_cast< PWORD * >( data + rIdx1 * blocks_per_row + 1 );
    PWORD *mb_tmp_data = reinterpret_cast< PWORD * >( data + rIdx2 * blocks_per_row + 1 );

    contingency_table case_cont, ctrl_cont;
    ResetContingencyTable( case_cont );
    ResetContingencyTable( ctrl_cont );

    register PWORD a_val, b_val, tmp_a, tmp_b, tmp_val;
    register PWORD a_aa, a_bb, b_aa, b_bb, a_ab, b_ab;
    PWORD mask = 0, uhalf, lhalf;

    const PWORD *case_ptr = reinterpret_cast< const PWORD * >( ccs.case_begin() );
    const PWORD *ctrl_ptr = reinterpret_cast< const PWORD * >( ccs.control_begin() );

    for( uint i = 1; i < blocks_per_row; i += BLOCKS_PER_PWORD ) {
        a_val = *ma_tmp_data;
        ++ma_tmp_data;
        b_val = *mb_tmp_data;
        ++mb_tmp_data;

        mask = *case_ptr++;
        if( mask ) {
            tmp_a = a_val & mask;
            tmp_b = b_val & mask;

            DecodeBitStrings2BitBlock( tmp_a, a_aa, a_ab, a_bb );
            DecodeBitStrings2BitBlock( tmp_b, b_aa, b_ab, b_bb );

            AddToContingency2BitBlock( case_cont.n0, tmp_val, a_aa, b_aa );
            AddToContingency2BitBlock( case_cont.n1, tmp_val, a_aa, b_ab );
            AddToContingency2BitBlock( case_cont.n2, tmp_val, a_aa, b_bb );
            AddToContingency2BitBlock( case_cont.n3, tmp_val, a_ab, b_aa );
            AddToContingency2BitBlock( case_cont.n4, tmp_val, a_ab, b_ab );
            AddToContingency2BitBlock( case_cont.n5, tmp_val, a_ab, b_bb );
            AddToContingency2BitBlock( case_cont.n6, tmp_val, a_bb, b_aa );
            AddToContingency2BitBlock( case_cont.n7, tmp_val, a_bb, b_ab );
            AddToContingency2BitBlock( case_cont.n8, tmp_val, a_bb, b_bb );
        }

        mask = *ctrl_ptr++;
        if( mask ) {
            tmp_a = a_val & mask;
            tmp_b = b_val & mask;

            DecodeBitStrings2BitBlock( tmp_a, a_aa, a_ab, a_bb );
            DecodeBitStrings2BitBlock( tmp_b, b_aa, b_ab, b_bb );

            AddToContingency2BitBlock( ctrl_cont.n0, tmp_val, a_aa, b_aa );
            AddToContingency2BitBlock( ctrl_cont.n1, tmp_val, a_aa, b_ab );
            AddToContingency2BitBlock( ctrl_cont.n2, tmp_val, a_aa, b_bb );
            AddToContingency2BitBlock( ctrl_cont.n3, tmp_val, a_ab, b_aa );
            AddToContingency2BitBlock( ctrl_cont.n4, tmp_val, a_ab, b_ab );
            AddToContingency2BitBlock( ctrl_cont.n5, tmp_val, a_ab, b_bb );
            AddToContingency2BitBlock( ctrl_cont.n6, tmp_val, a_bb, b_aa );
            AddToContingency2BitBlock( ctrl_cont.n7, tmp_val, a_bb, b_ab );
            AddToContingency2BitBlock( ctrl_cont.n8, tmp_val, a_bb, b_bb );
        }
    }
/*
    PWORD tail_end = ( max_column & TAIL_END );
    if( tail_end ) {
        a_val = *ma_tmp_data++;
        b_val = *mb_tmp_data++;

        tail_end = TAIL_END_MASK( tail_end );
        a_val &= tail_end;
        b_val &= tail_end;

        mask = *case_ptr++;
        mask &= tail_end;
        if( mask ) {
            tmp_a = a_val & mask;
            tmp_b = b_val & mask;

            DecodeBitStrings2BitBlock( tmp_a, a_aa, a_ab, a_bb );
            DecodeBitStrings2BitBlock( tmp_b, b_aa, b_ab, b_bb );

            AddToContingency2BitBlock( case_cont.n0, tmp_val, a_aa, b_aa );
            AddToContingency2BitBlock( case_cont.n1, tmp_val, a_aa, b_ab );
            AddToContingency2BitBlock( case_cont.n2, tmp_val, a_aa, b_bb );
            AddToContingency2BitBlock( case_cont.n3, tmp_val, a_ab, b_aa );
            AddToContingency2BitBlock( case_cont.n4, tmp_val, a_ab, b_ab );
            AddToContingency2BitBlock( case_cont.n5, tmp_val, a_ab, b_bb );
            AddToContingency2BitBlock( case_cont.n6, tmp_val, a_bb, b_aa );
            AddToContingency2BitBlock( case_cont.n7, tmp_val, a_bb, b_ab );
            AddToContingency2BitBlock( case_cont.n8, tmp_val, a_bb, b_bb );
        }

        mask = *ctrl_ptr++;
        mask &= tail_end;
        if( mask ) {
            tmp_a = a_val & mask;
            tmp_b = b_val & mask;

            DecodeBitStrings2BitBlock( tmp_a, a_aa, a_ab, a_bb );
            DecodeBitStrings2BitBlock( tmp_b, b_aa, b_ab, b_bb );

            AddToContingency2BitBlock( ctrl_cont.n0, tmp_val, a_aa, b_aa );
            AddToContingency2BitBlock( ctrl_cont.n1, tmp_val, a_aa, b_ab );
            AddToContingency2BitBlock( ctrl_cont.n2, tmp_val, a_aa, b_bb );
            AddToContingency2BitBlock( ctrl_cont.n3, tmp_val, a_ab, b_aa );
            AddToContingency2BitBlock( ctrl_cont.n4, tmp_val, a_ab, b_ab );
            AddToContingency2BitBlock( ctrl_cont.n5, tmp_val, a_ab, b_bb );
            AddToContingency2BitBlock( ctrl_cont.n6, tmp_val, a_bb, b_aa );
            AddToContingency2BitBlock( ctrl_cont.n7, tmp_val, a_bb, b_ab );
            AddToContingency2BitBlock( ctrl_cont.n8, tmp_val, a_bb, b_bb );
        }

    }
*/
    ccct.setMarkerAIndex( rIdx1 );
    ccct.setMarkerBIndex( rIdx2 );

    ccct.updateContingencyTables(case_cont, ctrl_cont);
}

//void CompressedGenotypeTable3::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct ) {
//    DataBlock *ma_tmp_data = data + rIdx1 * blocks_per_row + 1; // + 1 to skip header
//    DataBlock *mb_tmp_data = data + rIdx2 * blocks_per_row + 1; // + 1 to skip header
//
//    int tail_end = 0;
//    joint_genotypes val;
//    contingency_table case_cont, ctrl_cont;
//
//    ResetContingencyTable( case_cont );
//    ResetContingencyTable( ctrl_cont );
//
//    uint base_offset = 0xF0000;
//    ushort a_val, b_val;
//    const byte *case_ptr = ccs.case_begin(), *ctrl_ptr = ccs.control_begin();
//    byte _case, _ctrl;
//    uint case_skip_count = 0, ctrl_skip_count = 0;
//    ushort hi_offset, lo_offset;
//
//    for( uint i = 2; i < blocks_per_row; ++i, ++case_ptr, ++ctrl_ptr ) {
//        a_val = *ma_tmp_data++;
//        b_val = *mb_tmp_data++;
//
//        hi_offset = (( a_val & 0xFF00 ) | (( b_val & 0xFF00 ) >> 8 ) );
//        lo_offset = ((( a_val & 0x00FF ) << 8 ) | ( b_val & 0x00FF ) );
//
//        _case = *case_ptr;
//        _ctrl = *ctrl_ptr;
//
//        // compute case contingency
//        if( _case != 0 ) {
//            if(( base_offset = ( _case & 0xF0 ) ) != 0 ) {
//                case_skip_count += skip_count[ base_offset >> 4 ];
//                base_offset <<= 12;
//
//                val = contingency_lookup[ base_offset | hi_offset ];
//
//                ContingencyAddJointGenotype( case_cont, val );
//            }
//
//            if(( base_offset = ( _case & 0x0F ) ) != 0 ) {
//                case_skip_count += skip_count[ base_offset ];
//                base_offset <<= 16;
//
//                val = contingency_lookup[ base_offset | lo_offset ];
//                ContingencyAddJointGenotype( case_cont, val );
//            }
//        }
//
//        // compute control contingency
//        if( _ctrl != 0 ) {
//            if(( base_offset = ( _ctrl & 0xF0 ) ) != 0 ) {
//                ctrl_skip_count += skip_count[ base_offset >> 4 ];
//                base_offset <<= 12;
//
//                val = contingency_lookup[ base_offset | hi_offset ];
//
//                ContingencyAddJointGenotype( ctrl_cont, val );
//            }
//
//            if(( base_offset = ( _ctrl & 0x0F ) ) != 0 ) {
//
//                ctrl_skip_count += skip_count[ base_offset ];
//                base_offset <<= 16;
//
//                val = contingency_lookup[ base_offset | lo_offset ];
//
//                ContingencyAddJointGenotype( ctrl_cont, val );
//            }
//        }
//    }
//    // tail_end is the first index which is empty
//    tail_end = ( max_column & 7 );
//
//    a_val = *ma_tmp_data;
//    b_val = *mb_tmp_data;
//
//    _case = *case_ptr;
//    _ctrl = *ctrl_ptr;
//
//    while( tail_end > 0 ) {
//        hi_offset = (( a_val & 0xFF00 ) | (( b_val & 0xFF00 ) >> 8 ) );
//        if(( base_offset = ( _case & 0xF0 ) ) != 0 ) {
//            case_skip_count += skip_count[ base_offset >> 4 ];
//            base_offset <<= 12;
//            val = contingency_lookup[ base_offset | hi_offset ];
//
//            ContingencyAddJointGenotype( case_cont, val );
//        }
//
//        if(( base_offset = ( _ctrl & 0xF0 ) ) != 0 ) {
//            ctrl_skip_count += skip_count[ base_offset >> 4 ];
//            base_offset <<= 12;
//            val = contingency_lookup[ base_offset | hi_offset ];
//            ContingencyAddJointGenotype( ctrl_cont, val );
//        }
//
//        _case <<= 4;
//        _ctrl <<= 4;
//        a_val <<= 8;
//        b_val <<= 8;
//        tail_end -= 4;
//    }
//
//    case_cont.n9 -= case_skip_count;
//    ctrl_cont.n9 -= ctrl_skip_count;
//
//    ccct.setMarkerAIndex( rIdx1 );
//    ccct.setMarkerBIndex( rIdx2 );
//
//    ccct.setCaseContingency( case_cont );
//    ccct.setControlContingency( ctrl_cont );
//}

void CompressedGenotypeTable3::constructCountLookup( ) {
    genotype_counts counts;
    for( uint i = 0; i < 0x10000; ++i ) {
        counts.ui = 0;

        for( uint j = 0, k = i; j < 8; ++j, k >>= 2 ) {
            switch( k & 0x0003 ) {
            case 0x0003:
                ++counts.c3;
                break;
            case 0x0002:
                ++counts.c2;
                break;
            case 0x0001:
                ++counts.c1;
                break;
            default:
                ++counts.c0;
                break;
            }
        }

        count_lookup[i].ui = counts.ui;
    }
}

void CompressedGenotypeTable3::constructContingencyLookup() {
    joint_genotypes jg;

    byte ma, mb;
    uint base_offset, offset;
    int bit_mask, shift, idx;

    memset( contingency_lookup, 0, 16 * 256 * 256 * sizeof( joint_genotypes ) );

    for( uint k = 1, mask = 0x00; k < 16; ++k ) {
        mask = 0;
        mask |= (( k & 0x01 ) ? 0x03 : 0 );
        mask |= (( k & 0x02 ) ? 0x0C : 0 );
        mask |= (( k & 0x04 ) ? 0x30 : 0 );
        mask |= (( k & 0x08 ) ? 0xC0 : 0 );

        skip_count[ k ] = 0;
        skip_count[ k ] += ((( k & 0x01 ) ) ? 0 : 1 );
        skip_count[ k ] += ((( k & 0x02 ) ) ? 0 : 1 );
        skip_count[ k ] += ((( k & 0x04 ) ) ? 0 : 1 );
        skip_count[ k ] += ((( k & 0x08 ) ) ? 0 : 1 );

        for( uint i = 0; i < 256; ++i ) {
            base_offset = (( k << 16 ) | ( i << 8 ) );
            for( uint j = 0; j < 256; ++j ) {
                // clear data
                ResetJointGenotype( jg );
                for( bit_mask = 0xC0, shift = 6, idx = 0; shift >= 0; bit_mask >>= 2, shift -= 2, ++idx ) {
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

bool CompressedGenotypeTable3::isGenotypeHomozygous( ushort encoded_gt ) {
    // 0 == AA; 5 == CC; 10 == GG; 15 == TT
    return (( encoded_gt != 0xFFFF ) && ( encoded_gt == 0 || encoded_gt == 5 || encoded_gt == 10 || encoded_gt == 15 ) );
}

const char *CompressedGenotypeTable3::decodeGenotype( ushort encoded_gt ) {
//    return gt_lookup[encoded_gt];
    if( encoded_gt == 0xFFFF ) {
        return err_lookup;
    }
    return ( gt_lookup + encoded_gt * gt_size );
}

CompressedGenotypeTable3::~CompressedGenotypeTable3() {
    delete [] gt_lookup;

    for( uint i = 0; i < alphabet_size + 1; ++i ) {
        delete [] lookup[ i ];
    }
    delete [] lookup;

//    delete [] contingency_lookup;
//    delete [] skip_count;
}

}
}



