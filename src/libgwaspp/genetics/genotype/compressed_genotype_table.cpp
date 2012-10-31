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
#include "genetics/genotype/compressed_genotype_table.h"

namespace libgwaspp {
namespace genetics {

void CompressedGenotypeTable::initialize() {
    cout << "Initializing CompressedGenotypeTable ... " << endl;
    alphabet_size = ( int ) alphabet.length();
    possible_genotypes_size = pow( alphabet_size, ( double ) MAX_ALLELE_COUNT ) + 1;   // +1 for 0 -> "unknown" genotype

    cout << "Alphabet Size: " << alphabet_size << endl;
    cout << "Possible Genotype Size: " << possible_genotypes_size << endl;

    assert( possible_genotypes_size < ( 2 << BITS_PER_BLOCK ) );

    bits_per_data = 8;
    cout << "Bits per Data: " << bits_per_data << endl;

    assert(( BITS_PER_BLOCK % bits_per_data ) == 0 );   // assert that full data elements can be contained in a single block

    data_per_block = ( BITS_PER_BLOCK ) / bits_per_data;

    blocks_per_row = max_column / data_per_block + ((( max_column % data_per_block ) > 0 ) ? 1 : 0 );
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
    memset( transformations, ( byte )0x0F, 256 );

    int pattern = 1;
    for( string::const_iterator it = alphabet.begin(); it != alphabet.end(); ++it, pattern <<= 1 ) {
        transformations[( byte ) *it ] = pattern;
    }

    beg = new ushort[ possible_genotypes_size];
    memset( beg, 0xFF, possible_genotypes_size * sizeof( ushort ) );
    beg[ possible_genotypes_size - 1 ] = 0x00FF;
    end = beg + possible_genotypes_size;

    ushort val;
    for( uint i = 0, j = 0; j < alphabet_size; ++j ) {
        for( uint k = 0; k < alphabet_size; ++k ) {
            val = transformations[( byte ) alphabet[j] ];
            val <<= 4;
            val |= transformations[( byte ) alphabet[k] ];

            beg[i++] = val;
        }
    }
}

DataBlock CompressedGenotypeTable::operator()( int r, int c )  {
    // the multiplication step is a SLOW operation
    ulong offset = ( ulong ) r * blocks_per_row + (( uint ) c >> 1 );

//    ushort val = data[offset];
//    if( c & 0x01 ) {
//        val &= 0x00FF;
//    } else {
//        val &= 0xFF00;
//        val >>= 8;
//    }
    DataBlock db;
//    if( c & 0x01 ) {
//        db.us = data[offset].lo;
//    } else {
//        db.us = data[offset].hi;
//    }

    if( c & 0x01 ) {
        SetUshortAtDataBlock( db, GetLoByteAtDataBlock( data[offset] ));
    } else {
        SetUshortAtDataBlock( db, GetHiByteAtDataBlock( data[offset] ));
    }
    return db;
}

void CompressedGenotypeTable::addGenotype( int rIdx, int cIdx, const string &gt ) {
    assert( gt.length() == MAX_ALLELE_COUNT );

    int block_idx = cIdx / data_per_block;
    int block_offset = cIdx % data_per_block;

    DataBlock *tmp_data = data + rIdx * blocks_per_row + block_idx;

    ushort val;
    if( block_offset == 0 ) {
//        val = *tmp_data & 0x00FF;
        val = GetLoByteAtDataBlockPtr(tmp_data);

        val += (( ushort )transformations[( byte ) gt[0]] << 12 );
        val += (( ushort )transformations[( byte ) gt[0]] << 8 );
    } else {
//        val = *tmp_data & 0xFF00;
        val = GetHiByteAtDataBlockPtr(tmp_data);

        val += (( ushort )transformations[( byte ) gt[0]] << 4 );
        val += (( ushort )transformations[( byte ) gt[0]] );
    }

    SetUshortAtDataBlockPtr(tmp_data, val);
}

void CompressedGenotypeTable::addGenotypeRow( int rIdx, const char * p_begin, const char * p_end, char delim ) {
    assert(false);
}

void CompressedGenotypeTable::addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim ) {
#if MAX_ALLELE_COUNT == 2
    int cIdx = 0;

    ushort val = 0;
    DataBlock *tmp_data = data + rIdx * blocks_per_row;

    byte shift = 8, v0, v1;

    memset( tmp_data, ( byte ) 0xFF, bytes_per_row );

    do {
        if( shift > 8 ) {
//            *tmp_data = val;
//            tmp_data->us = val;
            SetUshortAtDataBlockPtr(tmp_data, val);
            ++tmp_data;
            val = 0;
            shift = 8;
        }
        if(( v0 = transformations[( byte ) *it++] ) == 0x0F ) {
            val |= ( ushort )( 0x00FF << shift );
            it++;
        } else if(( v1 = transformations[( byte ) *it++] ) == 0x0F ) {
            val |= ( ushort )( 0x00FF << shift );
        } else {
            val |= (( ushort )(( v0 << 4 ) | v1 ) << shift );
        }

        ++cIdx;
        shift -= 8;
    } while( it++ != it_end );

    SetUshortAtDataBlockPtr(tmp_data, val);
#else
#error Incomplete implementation of adding genotype by row
#endif
}

void CompressedGenotypeTable::selectMarker( uint rIdx ) {
    if( rIdx != current_dist_rIdx ) {
        resetSelectedMarker();

        current_dist_rIdx = rIdx;
        DataBlock *tmp_data = data + current_dist_rIdx * blocks_per_row;
        indices val;
        ushort code;

        memset( tmp_lookup, 0xFF, 0x10000 * sizeof( indices ) );
        tmp_lookup[ 0xFFFF ].us = 0xFFFF;  // both unknown

        ushort b0, b1;

        // stop at the last data block, as it may not be full
        for( uint i = 1; i < blocks_per_row; ++i, ++tmp_data ) {
            code = (ushort) *tmp_data;
            val = tmp_lookup[ code ];

            if( code == 0xFFFF ) {
                // two unknowns
                gt_dist.xx += 2;
            } else {
                if( val.us == 0xFFFF ) {
                    // new combination

                    // separate compacted genotypes
                    b0 = ((code >> 8) & 0x00FF);
                    b1 = ( code & 0x00FF );

                    if( b0 == 0x00FF ) {
                        val.hi = 0;
                    } else if( isGenotypeHomozygous( b0 ) ) {
                        if( gt_header.aa == 0 ) {
                            gt_header.aa = b0;
                            val.hi = 1;
                        } else if( gt_header.aa == b0 ) {
                            val.hi = 1;
                        } else if( gt_header.bb == 0 ) {
                            gt_header.bb = b0;
                            val.hi = 3;
                        } else if( gt_header.bb == b0 ) {
                            val.hi = 3;
                        } else {
                            assert( false );
                        }
                    } else {
                        if( gt_header.ab == 0 ) {
                            gt_header.ab = b0;
                        }
                        assert( gt_header.ab == b0 );
                        val.hi = 2;
                    }

                    if( b1 == 0x00FF ) {
                        val.lo = 0;
                    } else if( isGenotypeHomozygous( b1 ) ) {
                        if( gt_header.aa == 0 ) {
                            gt_header.aa = b1;
                            val.lo = 1;
                        } else if( gt_header.aa == b1 ) {
                            val.lo = 1;
                        } else if( gt_header.bb == 0 ) {
                            gt_header.bb = b1;
                            val.lo = 3;
                        } else if( gt_header.bb == b1 ) {
                            val.lo = 3;
                        } else {
                            assert( false );
                        }
                    } else {
                        if( gt_header.ab == 0 ) {
                            gt_header.ab = b1;
                        }
                        assert( gt_header.ab == b1 );
                        val.lo = 2;
                    }

                    tmp_lookup[ code ].us = val.us;

                    if( b0 != b1 ) {
                        tmp_lookup[(( b1 << 8 ) | b0 )].us = val.us;

//                        code = ( b0 << 8 ) | b0;
                        tmp_lookup[ ( b0 << 8 ) | b0 ].hi = val.hi;
                        tmp_lookup[ ( b0 << 8 ) | b0 ].lo = val.hi;

//                        code = ( b1 << 8 ) | b1;
                        tmp_lookup[ ( b1 << 8 ) | b1 ].hi = val.lo;
                        tmp_lookup[ ( b1 << 8 ) | b1 ].lo = val.lo;
                    }
                }

                ++gt_dist.freq[ val.hi ];
                ++gt_dist.freq[ val.lo ];
            }
        }

        // handle last data block independently of rest of blocks
        code = (ushort) *tmp_data;

        b0 = (( code >> 8 ) & 0x00FF );
        b1 = ( code & 0x00FF );

        if( b0 == 0x00FF ) {
            ++gt_dist.xx;
        } else if( isGenotypeHomozygous( b0 ) ) {
            if( gt_header.aa == 0 ) {
                gt_header.aa = b0;
                ++gt_dist.aa;
            } else if( gt_header.aa == b0 ) {
                ++gt_dist.aa;
            } else if( gt_header.bb == 0 ) {
                gt_header.bb = b0;
                ++gt_dist.bb;
            } else if( gt_header.bb == b0 ) {
                ++gt_dist.bb;
            } else {
                assert( false );
            }
        } else {
            if( gt_header.ab == 0 ) {
                gt_header.ab = b0;
            }
            assert( gt_header.ab == b0 );
            ++gt_dist.ab;
        }

        if( !( max_column & 1 ) ) {
            // if there are an even number of columns
            if( b1 == 0x00FF ) {
                ++gt_dist.xx;
            } else if( isGenotypeHomozygous( b1 ) ) {
                if( gt_header.aa == 0 ) {
                    gt_header.aa = b1;
                    ++gt_dist.aa;
                } else if( gt_header.aa == b1 ) {
                    ++gt_dist.aa;
                } else if( gt_header.bb == 0 ) {
                    gt_header.bb = b1;
                    ++gt_dist.bb;
                } else if( gt_header.bb == b1 ) {
                    ++gt_dist.bb;
                } else {
                    assert( false );
                }
            } else {
                if( gt_header.ab == 0 ) {
                    gt_header.ab = b1;
                }
                assert( gt_header.ab == b1 );
                ++gt_dist.ab;
            }
        }
        row_selected = true;
    }
}

void CompressedGenotypeTable::selectMarkerPair( uint rIdx1, uint rIdx2 ) {
    assert( false );
    if( ! ((maIdx == rIdx2 && mbIdx == rIdx1) || ( maIdx == rIdx1 && mbIdx == rIdx2 )) ) {
        resetSelectedMarkerPair();


    }
}

void CompressedGenotypeTable::getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist ) {
    assert( false );
}

void CompressedGenotypeTable::getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet &ccs, CaseControlGenotypeDistribution &ccgd ) {
    assert( false );
}

void CompressedGenotypeTable::getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct ) {
    assert( false );
}

void CompressedGenotypeTable::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct ) {
    assert( false );
}

ushort CompressedGenotypeTable::encodeGenotype( const string &gt ) {
    assert( gt.length() == MAX_ALLELE_COUNT );
    ushort val = transformations[( byte ) gt[0]];
    val <<= 4;
    val |= transformations[( byte ) gt[1]];

    return val;
}

bool CompressedGenotypeTable::isGenotypeHomozygous( ushort enc ) {
    return (( enc != 0xFFFF ) && (( enc & 0x000F ) == (( enc & 0x00F0 ) >> 4 ) ) );
}

const char *CompressedGenotypeTable::decodeGenotype( ushort encoded_gt ) {
    byte b0 = ( encoded_gt & 0x00F0 ), b1 = ( encoded_gt & 0x000F );
    if(( b1 == 0x000F ) || ( b0 == 0x00F0 ) || ( b0 == 0x0000 ) || ( b1 == 0x0000 ) ) {
        gt_lookup[0] = '0';
        gt_lookup[1] = '0';
    } else {
        int offset = -1;
        offset = (( b0 & 0x10 ) ? 0 : (( b0 & 0x20 ) ? 1 : (( b0 & 0x40 ) ? 2 : 3 ) ) );
        gt_lookup[0] = alphabet[offset];

        offset = (( b1 & 0x01 ) ? 0 : (( b1 & 0x02 ) ? 1 : (( b1 & 0x04 ) ? 2 : 3 ) ) );
        gt_lookup[1] = alphabet[offset];
    }
    return gt_lookup;
}

CompressedGenotypeTable::~CompressedGenotypeTable() {
    //dtor

    delete [] gt_lookup;
//    delete [] transformations;
}

}
}


