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
#include "genetics/genotype/compressed_genotype_table5.h"

namespace libgwaspp {
namespace genetics {

void CompressedGenotypeTable5::initialize() {
    cout << "Initializing CompressedGenotypeTable5 ... " << endl;
    alphabet_size = ( int ) alphabet.length();
    possible_genotypes_size = pow( alphabet_size, ( double ) MAX_ALLELE_COUNT ) + 1;   // +1 -> "unknown" genotype

    cout << "Alphabet Size: " << alphabet_size << endl;
    cout << "Possible Genotype Size: " << possible_genotypes_size << endl;

    assert( (int)possible_genotypes_size < ( 2 << BITS_PER_BLOCK ) );

    bits_per_data = 1;
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

    genotype_block_offset_ab = blocks_per_row;
    blocks_per_row *= 2;    // 2-bits per marker

    cout << "Genotype block offset: " << (int) genotype_block_offset_ab << endl;

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

void CompressedGenotypeTable5::addGenotype( int rIdx, int cIdx, const string &gt ) {
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

void CompressedGenotypeTable5::addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim ) {
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
    ushort val = 0, val_ab = 0;
    ushort shift = 1;

    ushort head_val = 0, geno_code = 0x0000;
    ushort clear_code = 0, head_shift = 0;

    ++tmp_data;         // skip header
    DataBlock *tmp_data_ab = tmp_data + genotype_block_offset_ab;
    do {
        if( shift == 0 ) {
            SetUshortAtDataBlockPtr( tmp_data, val );
            SetUshortAtDataBlockPtr( tmp_data_ab, val_ab );
            ++tmp_data;
            ++tmp_data_ab;
            val = 0;
            val_ab = 0;
            shift = 1;
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

            //val |= ( tmp_e << shift );
            switch( tmp_e ) {
            case 1:
                val |= shift;
                break;
            case 2:
                val_ab |= shift;
                break;
            case 3:
                val |= shift;
                val_ab |= shift;
                break;
            }
            tmp_e = 0;
        }

        shift <<=1;
    } while( ++it < it_end );

    SetUshortAtDataBlockPtr( tmp_data, val );
    SetUshortAtDataBlockPtr( tmp_data_ab, val_ab );
    SetUshortAtDataBlockPtr( tmp_header, head_val );
#else
#error "Incomplete implementation of adding genotype by row"
#endif
}

void CompressedGenotypeTable5::addGenotypeRow( int rIdx, const char *p_begin, const char *p_end, char delim ) {
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
    ushort val = 0, shift = 1;
    ushort val_ab = 0;

    ushort head_val = 0, geno_code = 0x0000;
    ushort clear_code = 0, head_shift = 0;

    ++tmp_data;
    const char *tmp_p = p_begin;
    DataBlock *tmp_data_ab = tmp_data + genotype_block_offset_ab;
    do {
        if( shift == 0 ) {
            SetUshortAtDataBlockPtr( tmp_data, val );
            SetUshortAtDataBlockPtr( tmp_data_ab, val_ab );
            ++tmp_data;
            ++tmp_data_ab;
            val = 0;
            val_ab = 0;
            shift = 1;
        }

        //cout << (char)*tmp_p;
        c1 = transformations[( byte )*tmp_p++];
        //cout << (char)*tmp_p;
        c2 = transformations[( byte )*tmp_p++];

        enc = GetUshortAtDataBlock( lookup[ c1 ][ c2 ] );
        //cout << " -> " << (int) enc;

        // only modify columns which are "known"
        if( enc != 0xFFFF ) {
            if(( tmp_e = enc_set[ enc ] ) == 0xFFFF ) {
                assert( geno_code < 0x7000 );

                HeaderStateMachine( geno_code, enc, head_shift, tmp_e, clear_code, head_val );

                enc_set[ enc ] = tmp_e;
            }

            assert( tmp_e != 0 );
            switch( tmp_e ) {
            case 1: // AA
                val |= shift;
                break;
            case 2: // AB
                val_ab |= shift;
                break;
            case 3: // BB
                val |= shift;
                val_ab |= shift;
                break;
            }
            /*cout << "(" << (int) tmp_e << ")";
            cout << hex;
            cout.width(5);
            cout << val;
            cout.width(5);
            cout << val_ab;
            cout << dec;*/
            tmp_e = 0;
        }
        //cout << endl;
        shift <<= 1;
    } while( ++tmp_p < p_end );

    SetUshortAtDataBlockPtr( tmp_data, val );
    SetUshortAtDataBlockPtr( tmp_data_ab, val_ab );
    // set header
    SetUshortAtDataBlockPtr( tmp_header, head_val );
#else
#error "Incomplete implementation of adding genotype by row"
#endif
}

ushort CompressedGenotypeTable5::encodeGenotype( const string &gt ) {
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

DataBlock CompressedGenotypeTable5::operator()( int r, int c ) {
    ulong row_offset = ( ulong ) r * blocks_per_row;
    ulong block_idx = (( uint ) c >> 4 );

    const ushort *header = data + row_offset;
    const ushort *tmp_data = header + block_idx + 1;

    ushort block_offset_mask = (1 << ( c & 0x0F ));

    DataBlock db;
    ushort v1 = ( *tmp_data );
    tmp_data += genotype_block_offset_ab;
    ushort v2 = (*tmp_data );

    v1 &= block_offset_mask;
    v2 &= block_offset_mask;

    if( v1 ) {
        if( v2 ) {  // BB
            SetUshortAtDataBlock( db, (*header & 0x000F) );
        } else {    // AA
            SetUshortAtDataBlock( db, ((*header & 0x0F00) >> 8) );
        }
    } else {
        if( v2 ) {  // AB
            SetUshortAtDataBlock( db, ((*header & 0x00F0) >> 4) );
        } else {
            SetUshortAtDataBlock( db, 0xFFFF );
        }
    }
    return db;
}


void CompressedGenotypeTable5::selectMarker( uint rIdx ) {
    assert(false);
}


void CompressedGenotypeTable5::selectMarkerPair( uint rIdx1, uint rIdx2 ) {
    assert(false);
}

void CompressedGenotypeTable5::selectCaseControl( CaseControlSet & ccs ) {
    if ( m_cases_controls == NULL ) {
        int tmp_data_per_block = ( BITS_PER_BLOCK ) / bits_per_data;

        // always pad the blocks per row by 1.
        // safe assumption that max_columns will not likely be evenly divisible by data_per_block 
        nCaseBlockCount = ccs.getCaseCount() / tmp_data_per_block + 1;
        nControlBlockCount = ccs.getControlCount() / tmp_data_per_block + 1;
        int block_bit_width = (sizeof(DataBlock) << 3);
        assert( (PROCESSOR_WORD_SIZE % block_bit_width) == 0 );

        // further pad the number of data blocks per row for processor word size efficiency
        int block_per_pword = (PROCESSOR_WORD_SIZE / block_bit_width);
        if( block_per_pword > 1) {
            if( nCaseBlockCount % block_per_pword ) {
                nCaseBlockCount += ( block_per_pword - (nCaseBlockCount % block_per_pword));
            }
            if( nControlBlockCount % block_per_pword ) {
                nControlBlockCount += ( block_per_pword -(nControlBlockCount % block_per_pword));
            }
        }
        nControlBlockOffset = 2 * nCaseBlockCount;
        nCaseControlBlockCount = 2 * (nCaseBlockCount + nControlBlockCount);

        nCaseControlSize = nCaseControlBlockCount * max_row;
        m_cases_controls = new DataBlock[ nCaseControlSize ]; 
    }

    nCaseCount = ccs.getCaseCount();
    nControlCount = ccs.getControlCount();
    nIndivids = nCaseCount + nControlCount;

    // clear case control buffer
    memset( m_cases_controls, 0, nCaseControlSize * sizeof(DataBlock) );

    // Mask out cases and controls
    const PWORD *case_ptr = reinterpret_cast< const PWORD * >( ccs.stream_case_begin() );
    const PWORD *ctrl_ptr = reinterpret_cast< const PWORD * >( ccs.stream_control_begin() );

    const PWORD *tmp_case, *tmp_ctrl;

    PWORD case_word, case_word_ab, ctrl_word, ctrl_word_ab;
    PWORD _aa, _ab;
    PWORD case_mask, ctrl_mask, mask;
    PWORD _case, _ctrl;

    // for every row
    for( uint i = 0; i < (uint)max_row; ++i ) {
        // locate the start of the data segment for each row 
        PWORD *_data = reinterpret_cast< PWORD * >( data + i * blocks_per_row + 1 );
        PWORD *_data_ab = reinterpret_cast< PWORD *>( data + i * blocks_per_row + 1 + genotype_block_offset_ab);

        PWORD *case_out = reinterpret_cast< PWORD * >( m_cases_controls + i * nCaseControlBlockCount );
        PWORD *case_out_ab = reinterpret_cast< PWORD * >(m_cases_controls + i * nCaseControlBlockCount + nCaseBlockCount );

        PWORD *ctrl_out = reinterpret_cast< PWORD * >( m_cases_controls + i * nCaseControlBlockCount + nControlBlockOffset );
        PWORD *ctrl_out_ab = reinterpret_cast< PWORD * >( m_cases_controls + i * nCaseControlBlockCount + nControlBlockOffset + nControlBlockCount);
        

        // reset the mask pointers
        tmp_case = case_ptr;
        tmp_ctrl = ctrl_ptr;

        case_word = 0;
        case_word_ab = 0;

        ctrl_word = 0;
        ctrl_word_ab = 0;

        case_mask = 1;
        ctrl_mask = 1;
        // for every data block
        for( uint j = 1; j < blocks_per_row; j += BLOCKS_PER_PWORD ) {
            // mask out all unnecessary data columns
            _aa = *_data++;
            _ab = *_data_ab++;

            _case = *tmp_case++;
            _ctrl = *tmp_ctrl++;

            mask = 1;
            while( mask != 0 ) {
                if( case_mask == 0 ) {
                    *case_out++ = case_word;
                    *case_out_ab++ = case_word_ab;

                    case_word = 0;
                    case_word_ab = 0;
                    case_mask = 1;
                }

                if( ctrl_mask == 0 ) {
                    *ctrl_out++ = ctrl_word;
                    *ctrl_out_ab++ = ctrl_word_ab;
                    ctrl_word = 0;
                    ctrl_word_ab = 0;
                    ctrl_mask = 1;
                }
                if( _case & mask ) {
                    if( _aa & mask ) {
                        case_word |= case_mask;
                        if( _ab & mask )
                            case_word_ab |= case_mask;
                    } else if( _ab & mask ) {
                        case_word_ab |= case_mask;
                    }
                    case_mask <<= 1;
                } else if( _ctrl & mask ) {
                    if (_aa & mask ) {
                        ctrl_word |= ctrl_mask;
                        if( _ab & mask )
                            ctrl_word_ab |= ctrl_mask;
                    } else if( _ab & mask ) {
                        ctrl_word_ab |= ctrl_mask;
                    }
                    ctrl_mask <<= 1;
                }

                mask <<= 1;
            }
            
        }

        if( case_mask != 1 ) {
            *case_out++ = case_word;
            *case_out_ab++ = case_word_ab;
        }
        if( ctrl_mask != 1 ) {
            *ctrl_out++ = ctrl_word;
            *ctrl_out_ab++ = ctrl_word_ab;
        }
    }
}

void CompressedGenotypeTable5::getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist ) {
    //ulong row_offset = rIdx * blocks_per_row + 1;
    PWORD *tmp_data = reinterpret_cast< PWORD * >(data + rIdx * blocks_per_row + 1);
    PWORD *tmp_data_ab = reinterpret_cast< PWORD * >(data +  rIdx * blocks_per_row + 1+ genotype_block_offset_ab);
    PWORD _aa, _ab;

    frequency_table joint_gt;
    ResetFrequencyTable(joint_gt);

    for( uint i = 1; i < genotype_block_offset_ab; i += BLOCKS_PER_PWORD ) {
        _aa = *tmp_data++;
        _ab = *tmp_data_ab++;

        /*cout << hex;
        cout.width(16);
        cout << _aa;
        cout << " -> " << dec << PopCount(_aa) << endl;*/
        //DecodeBitStreams2BitStream(_aa, _ab, _bb);
        //IncrementFrequencyValueStream(joint_gt, _aa, _ab, _bb);
        //
        joint_gt.bb += PopCount( _aa & _ab );
        joint_gt.ab += PopCount( _ab );
        joint_gt.aa += PopCount( _aa );
        joint_gt.xx += PopCount( ~( _aa | _ab ) );
    }

    joint_gt.aa -= joint_gt.bb;
    joint_gt.ab -= joint_gt.bb;

    dist.setDistribution( joint_gt );
}

void CompressedGenotypeTable5::getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet &ccs, CaseControlGenotypeDistribution &ccgd ) {
    //ulong row_offset = rIdx * blocks_per_row + 1;
    PWORD *tmp_data = reinterpret_cast< PWORD * >(data + rIdx * blocks_per_row + 1);
    PWORD *tmp_data_ab = reinterpret_cast< PWORD * >(data +  rIdx * blocks_per_row + 1+ genotype_block_offset_ab);

    register PWORD _aa, _ab;
    PWORD mask;

    frequency_table case_gt, ctrl_gt;
    ResetFrequencyTable( case_gt );
    ResetFrequencyTable( ctrl_gt );

    const PWORD *case_ptr = reinterpret_cast< const PWORD * >( ccs.stream_case_begin() );
    const PWORD *ctrl_ptr = reinterpret_cast< const PWORD * >( ccs.stream_control_begin() );

    // skip header block
    for( uint i = 1; i < genotype_block_offset_ab; i += BLOCKS_PER_PWORD ) {
        mask = *case_ptr++;
        _aa = mask & *tmp_data;
        _ab = mask & *tmp_data_ab;

        //DecodeBitStreams2BitStream( _aa, _ab, _bb );
        //IncrementFrequencyValueStream(case_gt, _aa, _ab, _bb);
        case_gt.bb += PopCount( _aa & _ab );
        case_gt.ab += PopCount( _ab );
        case_gt.aa += PopCount( _aa );

        mask = *ctrl_ptr++;
        _aa = mask & *tmp_data++;
        _ab = mask & *tmp_data_ab++;
        
        //DecodeBitStreams2BitStream( _aa, _ab, _bb );
        //IncrementFrequencyValueStream(ctrl_gt, _aa, _ab, _bb);
        ctrl_gt.bb += PopCount( _aa & _ab );
        ctrl_gt.ab += PopCount( _ab );
        ctrl_gt.aa += PopCount( _aa );
    }

    ctrl_gt.aa -= ctrl_gt.bb;
    ctrl_gt.ab -= ctrl_gt.bb;
    ctrl_gt.xx = ccs.getControlCount() - ctrl_gt.aa - ctrl_gt.ab - ctrl_gt.bb;

    case_gt.aa -= case_gt.bb;
    case_gt.ab -= case_gt.bb;
    case_gt.xx = ccs.getCaseCount() - case_gt.aa - case_gt.ab - case_gt.bb;

    ccgd.setCaseDistribution(case_gt);
    ccgd.setControlDistribution(ctrl_gt);
}

void CompressedGenotypeTable5::getCaseControlGenotypeDistribution( uint rIdx, CaseControlGenotypeDistribution &ccgd ) {
    //ulong offset = rIdx * nCaseControlBlockCount;
    PWORD *tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx * nCaseControlBlockCount);
    PWORD *tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx * nCaseControlBlockCount+ nCaseBlockCount);

    PWORD _aa, _ab;

    frequency_table case_gt, ctrl_gt;
    ResetFrequencyTable( case_gt );
    ResetFrequencyTable( ctrl_gt );

    for( uint i = 0; i < nCaseBlockCount; i += BLOCKS_PER_PWORD ) {
        _aa = *tmp_data++;
        _ab = *tmp_data_ab++;
        
        case_gt.bb += PopCount( _aa & _ab );
        case_gt.ab += PopCount( _ab );
        case_gt.aa += PopCount( _aa );
    }

    tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx * nCaseControlBlockCount+ nControlBlockOffset );
    tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx * nCaseControlBlockCount+ nControlBlockOffset + nControlBlockCount);

    for( uint i = 0; i < nControlBlockCount; i += BLOCKS_PER_PWORD ) {
        _aa = *tmp_data++;
        _ab = *tmp_data_ab++;

        ctrl_gt.bb += PopCount( _aa & _ab );
        ctrl_gt.ab += PopCount( _ab );
        ctrl_gt.aa += PopCount( _aa );
    }

    case_gt.ab -= case_gt.bb;
    case_gt.aa -= case_gt.bb;
    case_gt.xx = nCaseCount - case_gt.aa - case_gt.ab - case_gt.bb;

    ctrl_gt.ab -= ctrl_gt.bb;
    ctrl_gt.aa -= ctrl_gt.bb;
    ctrl_gt.xx = nControlCount - ctrl_gt.aa - ctrl_gt.ab - ctrl_gt.bb;

    ccgd.setCaseDistribution(case_gt);
    ccgd.setControlDistribution(ctrl_gt);
}

void CompressedGenotypeTable5::getCaseControlGenotypeDistribution( uint rIdx, CaseControlGenotypeDistribution &ccgd, marginal_information & m ) {
    //ulong offset = rIdx * nCaseControlBlockCount;
    PWORD *tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx * nCaseControlBlockCount);
    PWORD *tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx * nCaseControlBlockCount+ nCaseBlockCount);

    PWORD _aa, _ab;

    frequency_table case_gt, ctrl_gt;
    ResetFrequencyTable( case_gt );
    ResetFrequencyTable( ctrl_gt );

    for( uint i = 0; i < nCaseBlockCount; i += BLOCKS_PER_PWORD ) {
        _aa = *tmp_data++;
        _ab = *tmp_data_ab++;
        
        case_gt.bb += PopCount( _aa & _ab );
        case_gt.ab += PopCount( _ab );
        case_gt.aa += PopCount( _aa );
    }

    tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx * nCaseControlBlockCount+ nControlBlockOffset );
    tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx * nCaseControlBlockCount+ nControlBlockOffset + nControlBlockCount);

    for( uint i = 0; i < nControlBlockCount; i += BLOCKS_PER_PWORD ) {
        _aa = *tmp_data++;
        _ab = *tmp_data_ab++;

        ctrl_gt.bb += PopCount( _aa & _ab );
        ctrl_gt.ab += PopCount( _ab );
        ctrl_gt.aa += PopCount( _aa );
    }

    case_gt.ab -= case_gt.bb;
    case_gt.aa -= case_gt.bb;
    case_gt.xx = nCaseCount - case_gt.aa - case_gt.ab - case_gt.bb;

    ctrl_gt.ab -= ctrl_gt.bb;
    ctrl_gt.aa -= ctrl_gt.bb;
    ctrl_gt.xx = nControlCount - ctrl_gt.aa - ctrl_gt.ab - ctrl_gt.bb;

    ccgd.setCaseDistribution(case_gt);
    ccgd.setControlDistribution(ctrl_gt);

    computeMarginalInformation( case_gt, ctrl_gt, (double)nIndivids, m );
}

void CompressedGenotypeTable5::getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct ) {
    //ulong ma_offset = rIdx1 * blocks_per_row + 1;
    //ulong mb_offset = rIdx2 * blocks_per_row + 1;
    PWORD *ma_tmp_data = reinterpret_cast< PWORD * >( data +  rIdx1 * blocks_per_row + 1);
    PWORD *ma_tmp_data_ab = reinterpret_cast< PWORD * >( data +  rIdx1 * blocks_per_row + 1+ genotype_block_offset_ab );

    PWORD *mb_tmp_data = reinterpret_cast< PWORD * >( data +  rIdx2 * blocks_per_row + 1);
    PWORD *mb_tmp_data_ab = reinterpret_cast< PWORD * >( data +  rIdx2 * blocks_per_row + 1+ genotype_block_offset_ab );

    ct.setMarkerAIndex( rIdx1 );
    ct.setMarkerBIndex( rIdx2 );

    CONTIN_TABLE_T contingency;
    ResetContingencyTable( contingency );

    PWORD a_aa, a_bb, b_aa, b_bb, a_ab, b_ab, a_xx, b_xx;

    for( uint i = 1; i < genotype_block_offset_ab; i += BLOCKS_PER_PWORD ) {
        a_aa = *ma_tmp_data++;
        a_ab = *ma_tmp_data_ab++;
        
        b_aa = *mb_tmp_data++;
        b_ab = *mb_tmp_data_ab++;

        DecodeBitStreams2BitStream( a_aa, a_ab, a_bb, a_xx );
        DecodeBitStreams2BitStream( b_aa, b_ab, b_bb, b_xx );

        AddToContingencyStream( contingency.AA_BB, a_aa, b_aa );
        AddToContingencyStream( contingency.AA_Bb, a_aa, b_ab );
        AddToContingencyStream( contingency.AA_bb, a_aa, b_bb );
        AddToContingencyStream( contingency.Aa_BB, a_ab, b_aa );
        AddToContingencyStream( contingency.Aa_Bb, a_ab, b_ab );
        AddToContingencyStream( contingency.Aa_bb, a_ab, b_bb );
        AddToContingencyStream( contingency.aa_BB, a_bb, b_aa );
        AddToContingencyStream( contingency.aa_Bb, a_bb, b_ab );
        AddToContingencyStream( contingency.aa_bb, a_bb, b_bb );

        if( b_xx ) {
            AddToContingencyStream( contingency.AA_xx, a_aa, b_xx );
            AddToContingencyStream( contingency.Aa_xx, a_ab, b_xx );
            AddToContingencyStream( contingency.aa_xx, a_bb, b_xx );
        }
        if( a_xx ) {
            AddToContingencyStream( contingency.xx_BB, a_xx, b_aa );
            AddToContingencyStream( contingency.xx_Bb, a_xx, b_ab );
            AddToContingencyStream( contingency.xx_bb, a_xx, b_bb );
        }
        if( a_xx || b_xx )
            AddToContingencyStream( contingency.xx_xx, a_xx, b_xx );
    }
    ct.setContingency( contingency );
}

void CompressedGenotypeTable5::getContingencyTable( uint rIdx1, uint rIdx2, ushort *column_set, ContingencyTable &ct ) {
    assert(false);
}

void CompressedGenotypeTable5::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct ) {
    //ulong ma_offset = rIdx1 * blocks_per_row + 1;
    //ulong mb_offset = rIdx2 * blocks_per_row + 1;
    PWORD *ma_tmp_data = reinterpret_cast< PWORD * >( data +  rIdx1 * blocks_per_row + 1);
    PWORD *ma_tmp_data_ab = reinterpret_cast< PWORD * >( data +  rIdx1 * blocks_per_row + 1+ genotype_block_offset_ab );

    PWORD *mb_tmp_data = reinterpret_cast< PWORD * >( data +  rIdx2 * blocks_per_row + 1);
    PWORD *mb_tmp_data_ab = reinterpret_cast< PWORD * >( data +  rIdx2 * blocks_per_row + 1+ genotype_block_offset_ab );

    CONTIN_TABLE_T case_cont, ctrl_cont;
    ResetContingencyTable( case_cont );
    ResetContingencyTable( ctrl_cont );

    register PWORD a_aa, a_bb, b_aa, b_bb, a_ab, b_ab, b_xx, a_xx;
    PWORD mask = 0;

    const PWORD *case_ptr = reinterpret_cast< const PWORD * >( ccs.stream_case_begin() );
    const PWORD *ctrl_ptr = reinterpret_cast< const PWORD * >( ccs.stream_control_begin() );

    for( uint i = 1; i < genotype_block_offset_ab; i += BLOCKS_PER_PWORD ) {
        mask = *case_ptr++;

        a_aa = mask & *ma_tmp_data;
        a_ab = mask & *ma_tmp_data_ab;
        DecodeBitStreams2BitStream( a_aa, a_ab, a_bb, a_xx);

        b_aa = mask & *mb_tmp_data;
        b_ab = mask & *mb_tmp_data_ab;
        DecodeBitStreams2BitStream( b_aa, b_ab, b_bb, b_xx);

        AddToContingencyStream( case_cont.AA_BB, a_aa, b_aa );
        AddToContingencyStream( case_cont.AA_Bb, a_aa, b_ab );
        AddToContingencyStream( case_cont.AA_bb, a_aa, b_bb );
        AddToContingencyStream( case_cont.Aa_BB, a_ab, b_aa );
        AddToContingencyStream( case_cont.Aa_Bb, a_ab, b_ab );
        AddToContingencyStream( case_cont.Aa_bb, a_ab, b_bb );
        AddToContingencyStream( case_cont.aa_BB, a_bb, b_aa );
        AddToContingencyStream( case_cont.aa_Bb, a_bb, b_ab );
        AddToContingencyStream( case_cont.aa_bb, a_bb, b_bb );

        if( b_xx ) {
            AddToContingencyStream( case_cont.AA_xx, a_aa, b_xx );
            AddToContingencyStream( case_cont.Aa_xx, a_ab, b_xx );
            AddToContingencyStream( case_cont.aa_xx, a_bb, b_xx );
        }
        if( a_xx ) {
            AddToContingencyStream( case_cont.xx_BB, a_xx, b_aa );
            AddToContingencyStream( case_cont.xx_Bb, a_xx, b_ab );
            AddToContingencyStream( case_cont.xx_bb, a_xx, b_bb );
        }
        if( a_xx || b_xx )
            AddToContingencyStream( case_cont.xx_xx, a_xx, b_xx );

        mask = *ctrl_ptr++;
        a_aa = mask & *ma_tmp_data++;
        a_ab = mask & *ma_tmp_data_ab++;
        DecodeBitStreams2BitStream( a_aa, a_ab, a_bb, a_xx );

        b_aa = mask & *mb_tmp_data++;
        b_ab = mask & *mb_tmp_data_ab++;
        DecodeBitStreams2BitStream( b_aa, b_ab, b_bb, b_xx );

        AddToContingencyStream( ctrl_cont.AA_BB, a_aa, b_aa );
        AddToContingencyStream( ctrl_cont.AA_Bb, a_aa, b_ab );
        AddToContingencyStream( ctrl_cont.AA_bb, a_aa, b_bb );
        AddToContingencyStream( ctrl_cont.Aa_BB, a_ab, b_aa );
        AddToContingencyStream( ctrl_cont.Aa_Bb, a_ab, b_ab );
        AddToContingencyStream( ctrl_cont.Aa_bb, a_ab, b_bb );
        AddToContingencyStream( ctrl_cont.aa_BB, a_bb, b_aa );
        AddToContingencyStream( ctrl_cont.aa_Bb, a_bb, b_ab );
        AddToContingencyStream( ctrl_cont.aa_bb, a_bb, b_bb );

        if( b_xx ) {
            AddToContingencyStream( ctrl_cont.AA_xx, a_aa, b_xx );
            AddToContingencyStream( ctrl_cont.Aa_xx, a_ab, b_xx );
            AddToContingencyStream( ctrl_cont.aa_xx, a_bb, b_xx );
        }
        if( a_xx ) {
            AddToContingencyStream( ctrl_cont.xx_BB, a_xx, b_aa );
            AddToContingencyStream( ctrl_cont.xx_Bb, a_xx, b_ab );
            AddToContingencyStream( ctrl_cont.xx_bb, a_xx, b_bb );
        }
        if( a_xx || b_xx )
            AddToContingencyStream( ctrl_cont.xx_xx, a_xx, b_xx );
    }
    ccct.setMarkerAIndex( rIdx1 );
    ccct.setMarkerBIndex( rIdx2 );

    ccct.updateContingencyTables(case_cont, ctrl_cont);
}
void CompressedGenotypeTable5::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlContingencyTable &ccct ) {
    //ulong ma_offset = rIdx1 * nCaseControlBlockCount;
    //ulong mb_offset = rIdx2 * nCaseControlBlockCount;
    PWORD *ma_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx1 * nCaseControlBlockCount);
    PWORD *ma_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx1 * nCaseControlBlockCount+ nCaseBlockCount );

    PWORD *mb_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx2 * nCaseControlBlockCount);
    PWORD *mb_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx2 * nCaseControlBlockCount+ nCaseBlockCount );

    CONTIN_TABLE_T case_cont, ctrl_cont;
    ResetContingencyTable( case_cont );
    ResetContingencyTable( ctrl_cont );

    PWORD a_aa, a_bb, b_aa, b_bb, a_ab, b_ab, a_xx, b_xx;

    for( uint i = 0; i < nCaseBlockCount; i += BLOCKS_PER_PWORD ) {
        a_aa = *ma_tmp_data++;
        a_ab = *ma_tmp_data_ab++;

        b_aa = *mb_tmp_data++;
        b_ab = *mb_tmp_data_ab++;

        DecodeBitStreams2BitStream( a_aa, a_ab, a_bb, a_xx );
        DecodeBitStreams2BitStream( b_aa, b_ab, b_bb, b_xx );

        AddToContingencyStream( case_cont.AA_BB, a_aa, b_aa );
        AddToContingencyStream( case_cont.AA_Bb, a_aa, b_ab );
        AddToContingencyStream( case_cont.AA_bb, a_aa, b_bb );
        AddToContingencyStream( case_cont.Aa_BB, a_ab, b_aa );
        AddToContingencyStream( case_cont.Aa_Bb, a_ab, b_ab );
        AddToContingencyStream( case_cont.Aa_bb, a_ab, b_bb );
        AddToContingencyStream( case_cont.aa_BB, a_bb, b_aa );
        AddToContingencyStream( case_cont.aa_Bb, a_bb, b_ab );
        AddToContingencyStream( case_cont.aa_bb, a_bb, b_bb );
        if( b_xx ) {
            AddToContingencyStream( case_cont.AA_xx, a_aa, b_xx );
            AddToContingencyStream( case_cont.Aa_xx, a_ab, b_xx );
            AddToContingencyStream( case_cont.aa_xx, a_bb, b_xx );
        }
        if( a_xx ) {
            AddToContingencyStream( case_cont.xx_BB, a_xx, b_aa );
            AddToContingencyStream( case_cont.xx_Bb, a_xx, b_ab );
            AddToContingencyStream( case_cont.xx_bb, a_xx, b_bb );
        }
        if( a_xx || b_xx )
            AddToContingencyStream( case_cont.xx_xx, a_xx, b_xx );
    }

    ma_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx1 * nCaseControlBlockCount+ nControlBlockOffset );
    ma_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx1 * nCaseControlBlockCount+ nControlBlockOffset + nControlBlockCount );

    mb_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx2 * nCaseControlBlockCount+ nControlBlockOffset );
    mb_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx2 * nCaseControlBlockCount+ nControlBlockOffset + nControlBlockCount );

    for( uint i = 0; i < nControlBlockCount; i += BLOCKS_PER_PWORD ) {
        a_aa = *ma_tmp_data++;
        a_ab = *ma_tmp_data_ab++;

        b_aa = *mb_tmp_data++;
        b_ab = *mb_tmp_data_ab++;

        DecodeBitStreams2BitStream( a_aa, a_ab, a_bb, a_xx );
        DecodeBitStreams2BitStream( b_aa, b_ab, b_bb, b_xx );

        AddToContingencyStream( ctrl_cont.AA_BB, a_aa, b_aa );
        AddToContingencyStream( ctrl_cont.AA_Bb, a_aa, b_ab );
        AddToContingencyStream( ctrl_cont.AA_bb, a_aa, b_bb );
        AddToContingencyStream( ctrl_cont.Aa_BB, a_ab, b_aa );
        AddToContingencyStream( ctrl_cont.Aa_Bb, a_ab, b_ab );
        AddToContingencyStream( ctrl_cont.Aa_bb, a_ab, b_bb );
        AddToContingencyStream( ctrl_cont.aa_BB, a_bb, b_aa );
        AddToContingencyStream( ctrl_cont.aa_Bb, a_bb, b_ab );
        AddToContingencyStream( ctrl_cont.aa_bb, a_bb, b_bb );
        if( b_xx ) {
            AddToContingencyStream( ctrl_cont.AA_xx, a_aa, b_xx );
            AddToContingencyStream( ctrl_cont.Aa_xx, a_ab, b_xx );
            AddToContingencyStream( ctrl_cont.aa_xx, a_bb, b_xx );
        }
        if( a_xx ) {
            AddToContingencyStream( ctrl_cont.xx_BB, a_xx, b_aa );
            AddToContingencyStream( ctrl_cont.xx_Bb, a_xx, b_ab );
            AddToContingencyStream( ctrl_cont.xx_bb, a_xx, b_bb );
        }
        if( a_xx || b_xx )
            AddToContingencyStream( ctrl_cont.xx_xx, a_xx, b_xx );
    }

    ccct.setMarkerAIndex( rIdx1 );
    ccct.setMarkerBIndex( rIdx2 );

    ccct.updateContingencyTables(case_cont, ctrl_cont);
}

void CompressedGenotypeTable5::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, const marginal_information &m1, const marginal_information &m2, CaseControlContingencyTable & ccct ) {
    PWORD *ma_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx1 * nCaseControlBlockCount);
    PWORD *ma_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls + rIdx1 * nCaseControlBlockCount + nCaseBlockCount );

    PWORD *mb_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls +  rIdx2 * nCaseControlBlockCount);
    PWORD *mb_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls + rIdx2 * nCaseControlBlockCount+ nCaseBlockCount );

    CONTIN_TABLE_T case_cont, ctrl_cont;
    ResetContingencyTable( case_cont );
    ResetContingencyTable( ctrl_cont );

    if( m1.cases.xx + m1.controls.xx + m2.cases.xx + m2.controls.xx ) {
        register PWORD a_aa, b_aa, a_bb, b_bb, a_ab, b_ab;
        for( uint i = 0; i < nCaseBlockCount; i += BLOCKS_PER_PWORD ) {
            a_aa = *ma_tmp_data++;
            a_ab = *ma_tmp_data_ab++;

            b_aa = *mb_tmp_data++;
            b_ab = *mb_tmp_data_ab++;

            DecodeBitStreams2BitStream( a_aa, a_ab, a_bb );
            DecodeBitStreams2BitStream( b_aa, b_ab, b_bb );

            AddToContingencyStream( case_cont.AA_BB, a_aa, b_aa );
            AddToContingencyStream( case_cont.AA_Bb, a_aa, b_ab );
            AddToContingencyStream( case_cont.AA_bb, a_aa, b_bb );
            AddToContingencyStream( case_cont.Aa_BB, a_ab, b_aa );
            AddToContingencyStream( case_cont.Aa_Bb, a_ab, b_ab );
            AddToContingencyStream( case_cont.Aa_bb, a_ab, b_bb );
            AddToContingencyStream( case_cont.aa_BB, a_bb, b_aa );
            AddToContingencyStream( case_cont.aa_Bb, a_bb, b_ab );
            AddToContingencyStream( case_cont.aa_bb, a_bb, b_bb );
        }
        case_cont.AA_xx = m1.cases.aa - case_cont.AA_BB - case_cont.AA_bb - case_cont.AA_Bb;
        case_cont.Aa_xx = m1.cases.ab - case_cont.Aa_BB - case_cont.Aa_bb - case_cont.Aa_Bb;
        case_cont.aa_xx = m1.cases.bb - case_cont.aa_bb - case_cont.aa_BB - case_cont.aa_Bb;

        case_cont.xx_BB = m2.cases.aa - case_cont.AA_BB - case_cont.aa_BB - case_cont.Aa_BB;
        case_cont.xx_Bb = m2.cases.ab - case_cont.AA_Bb - case_cont.aa_Bb - case_cont.Aa_Bb;
        case_cont.xx_bb = m2.cases.bb - case_cont.AA_bb - case_cont.aa_bb - case_cont.Aa_bb;

        case_cont.xx_xx = m2.cases.xx - case_cont.AA_xx - case_cont.Aa_xx - case_cont.aa_xx;

        ma_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls + rIdx1 * nCaseControlBlockCount + nControlBlockOffset );
        ma_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls + rIdx1 * nCaseControlBlockCount+ nControlBlockOffset + nControlBlockCount );

        mb_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls + rIdx2 * nCaseControlBlockCount + nControlBlockOffset );
        mb_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls + rIdx2 * nCaseControlBlockCount + nControlBlockOffset + nControlBlockCount );

        for( uint i = 0; i < nControlBlockCount; i += BLOCKS_PER_PWORD ) {
            a_aa = *ma_tmp_data++;
            a_ab = *ma_tmp_data_ab++;

            b_aa = *mb_tmp_data++;
            b_ab = *mb_tmp_data_ab++;

            DecodeBitStreams2BitStream( a_aa, a_ab, a_bb );
            DecodeBitStreams2BitStream( b_aa, b_ab, b_bb );

            AddToContingencyStream( ctrl_cont.AA_BB, a_aa, b_aa );
            AddToContingencyStream( ctrl_cont.AA_Bb, a_aa, b_ab );
            AddToContingencyStream( ctrl_cont.AA_bb, a_aa, b_bb );
            AddToContingencyStream( ctrl_cont.Aa_BB, a_ab, b_aa );
            AddToContingencyStream( ctrl_cont.Aa_Bb, a_ab, b_ab );
            AddToContingencyStream( ctrl_cont.Aa_bb, a_ab, b_bb );
            AddToContingencyStream( ctrl_cont.aa_BB, a_bb, b_aa );
            AddToContingencyStream( ctrl_cont.aa_Bb, a_bb, b_ab );
            AddToContingencyStream( ctrl_cont.aa_bb, a_bb, b_bb );
        }

        ctrl_cont.AA_xx = m1.controls.aa - ctrl_cont.AA_BB - ctrl_cont.AA_bb - ctrl_cont.AA_Bb;
        ctrl_cont.Aa_xx = m1.controls.ab - ctrl_cont.Aa_BB - ctrl_cont.Aa_bb - ctrl_cont.Aa_Bb;
        ctrl_cont.aa_xx = m1.controls.bb - ctrl_cont.aa_bb - ctrl_cont.aa_BB - ctrl_cont.aa_Bb;

        ctrl_cont.xx_BB = m2.controls.aa - ctrl_cont.AA_BB - ctrl_cont.aa_BB - ctrl_cont.Aa_BB;
        ctrl_cont.xx_Bb = m2.controls.ab - ctrl_cont.AA_Bb - ctrl_cont.aa_Bb - ctrl_cont.Aa_Bb;
        ctrl_cont.xx_bb = m2.controls.bb - ctrl_cont.AA_bb - ctrl_cont.aa_bb - ctrl_cont.Aa_bb;

        ctrl_cont.xx_xx = m2.controls.xx - ctrl_cont.AA_xx - ctrl_cont.Aa_xx - ctrl_cont.aa_xx;

    } else {
        // shortcut
        // since the margins for both markers are known
        // if there are no individuals with missing genotypes
        // then only the AA_BB, AA_bb, and Aa_
        //
        register PWORD a_aa, b_aa, a_bb, b_bb;
        register PWORD count_AB = 0, count_Ab = 0, count_aB = 0, count_ab = 0;
        for( uint i = 0; i < nCaseBlockCount; i += BLOCKS_PER_PWORD ) {
            a_aa = *ma_tmp_data++;
            a_bb = *ma_tmp_data_ab++;

            b_aa = *mb_tmp_data++;
            b_bb = *mb_tmp_data_ab++;

            a_bb &= a_aa; // bb
            a_aa ^= a_bb;       // aa
            b_bb &= b_aa; // bb
            b_aa ^= b_bb;       // aa

            AddToContingencyStream( count_AB, a_aa, b_aa );
            AddToContingencyStream( count_Ab, a_aa, b_bb );
            AddToContingencyStream( count_aB, a_bb, b_aa );
            AddToContingencyStream( count_ab, a_bb, b_bb );
        }

        case_cont.AA_BB = count_AB;
        case_cont.AA_bb = count_Ab;
        case_cont.aa_BB = count_aB;
        case_cont.aa_bb = count_ab;

        case_cont.AA_Bb = m1.cases.aa - case_cont.AA_BB - case_cont.AA_bb;
        case_cont.aa_Bb = m1.cases.bb - case_cont.aa_bb - case_cont.aa_BB;

        case_cont.Aa_BB = m2.cases.aa - case_cont.AA_BB - case_cont.aa_BB;
        case_cont.Aa_bb = m2.cases.bb - case_cont.AA_bb - case_cont.aa_bb;

        case_cont.Aa_Bb = m2.cases.ab - case_cont.AA_Bb - case_cont.aa_Bb;

        ma_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls + rIdx1 * nCaseControlBlockCount + nControlBlockOffset );
        ma_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls + rIdx1 * nCaseControlBlockCount+ nControlBlockOffset + nControlBlockCount );

        mb_tmp_data = reinterpret_cast< PWORD * >( m_cases_controls + rIdx2 * nCaseControlBlockCount + nControlBlockOffset );
        mb_tmp_data_ab = reinterpret_cast< PWORD * >( m_cases_controls + rIdx2 * nCaseControlBlockCount + nControlBlockOffset + nControlBlockCount );

        count_AB = 0; count_Ab = 0; count_aB = 0; count_ab = 0;
        for( uint i = 0; i < nControlBlockCount; i += BLOCKS_PER_PWORD ) {
            a_aa = *ma_tmp_data++;
            a_bb = *ma_tmp_data_ab++;

            b_aa = *mb_tmp_data++;
            b_bb = *mb_tmp_data_ab++;

            a_bb &= a_aa;
            a_aa ^= a_bb;
            b_bb &= b_aa;
            b_aa ^= b_bb;

            AddToContingencyStream( count_AB, a_aa, b_aa );
            AddToContingencyStream( count_Ab, a_aa, b_bb );
            AddToContingencyStream( count_aB, a_bb, b_aa );
            AddToContingencyStream( count_ab, a_bb, b_bb );
        }
        ctrl_cont.AA_BB = count_AB;
        ctrl_cont.AA_bb = count_Ab;
        ctrl_cont.aa_BB = count_aB;
        ctrl_cont.aa_bb = count_ab;

        ctrl_cont.AA_Bb = m1.controls.aa - ctrl_cont.AA_BB - ctrl_cont.AA_bb;
        ctrl_cont.aa_Bb = m1.controls.bb - ctrl_cont.aa_bb - ctrl_cont.aa_BB;

        ctrl_cont.Aa_BB = m2.controls.aa - ctrl_cont.AA_BB - ctrl_cont.aa_BB;
        ctrl_cont.Aa_bb = m2.controls.bb - ctrl_cont.AA_bb - ctrl_cont.aa_bb;

        ctrl_cont.Aa_Bb = m2.controls.ab - ctrl_cont.AA_Bb - ctrl_cont.aa_Bb;
    }
    
    ccct.setMarkerAIndex( rIdx1 );
    ccct.setMarkerBIndex( rIdx2 );

    ccct.updateContingencyTables(case_cont, ctrl_cont);
}

void CompressedGenotypeTable5::constructCountLookup( ) {
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

void CompressedGenotypeTable5::constructContingencyLookup() {
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

bool CompressedGenotypeTable5::isGenotypeHomozygous( ushort encoded_gt ) {
    // 0 == AA; 5 == CC; 10 == GG; 15 == TT
    return (( encoded_gt != 0xFFFF ) && ( encoded_gt == 0 || encoded_gt == 5 || encoded_gt == 10 || encoded_gt == 15 ) );
}

const char *CompressedGenotypeTable5::decodeGenotype( ushort encoded_gt ) {
//    return gt_lookup[encoded_gt];
    if( encoded_gt == 0xFFFF ) {
        return err_lookup;
    }
    return ( gt_lookup + encoded_gt * gt_size );
}

CompressedGenotypeTable5::~CompressedGenotypeTable5() {
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



