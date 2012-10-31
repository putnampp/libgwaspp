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
#include "genetics/genotype/basic_genotype_table.h"

namespace libgwaspp {
namespace genetics {

void BasicGenotypeTable::initialize() {
    cout << "Initializing BasicGenotypeTable ... " << endl;
    alphabet_size = ( int ) alphabet.length();
    possible_genotypes_size = pow( alphabet_size, ( double ) MAX_ALLELE_COUNT ) + 1;   // +1 for 0 -> "unknown" genotype

    cout << "Alphabet Size: " << alphabet_size << endl;
    cout << "Possible Genotype Size: " << possible_genotypes_size << endl;

    assert( possible_genotypes_size < ( 2 << BITS_PER_BLOCK ) );

    bits_per_data = 16;
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

    memset( transformations, 0xFF, 256 );

    for( string::const_iterator it = alphabet.begin(); it != alphabet.end(); it++ ) {
        transformations[( byte ) *it ] = *it;
    }

    beg = new ushort[ possible_genotypes_size ];
    memset( beg, 0xFF, possible_genotypes_size * sizeof( ushort ) );
    end = beg + possible_genotypes_size;

    ushort val;
    for( uint i = 0, j = 0; j < alphabet_size; ++j ) {
        for( uint k = 0; k < alphabet_size; ++k ) {
            val = ( byte ) alphabet[j];
            val <<= 8;
            val |= ( byte ) alphabet[k];

            beg[i++] = val;
        }
    }
}

DataBlock BasicGenotypeTable::operator()( int r, int c ) {
    ulong offset = ( ulong )r * blocks_per_row + c;

    return data[offset];
}

bool BasicGenotypeTable::isGenotypeHomozygous( ushort enc ) {
    return (( enc != 0xFFFF ) && ((( enc & 0xFF00 ) >> 8 ) == ( enc & 0x00FF ) ) );
}

void BasicGenotypeTable::addGenotype( int rIdx, int cIdx, const string &gt ) {
    assert( gt.length() == MAX_ALLELE_COUNT );
    ushort val = 0;

    ushort v0, v1;

    if((( v0 = transformations[( byte ) gt[0] ] ) == 0xFF ) || (( v1 = transformations[( byte ) gt[1] ] ) == 0xFF ) ) {
        val = 0xFFFF;
    } else {
        val = ( v0 << 8 ) | v1;
    }

    ulong offset = ( ulong )rIdx * blocks_per_row + ( uint )cIdx;

    SetUshortAtDataBlock( data[ offset ], val );
}

void BasicGenotypeTable::addGenotypeRow( int rIdx, const char * p_begin, const char * p_end, char delim ) {
    assert(false);
}

void BasicGenotypeTable::addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim ) {
#if MAX_ALLELE_COUNT == 2
    int cIdx = 0;

    ushort val, v0, v1;
    DataBlock *tmp_data = data + rIdx * blocks_per_row;

    memset( tmp_data, 0, bytes_per_row );

    do {
        if(( v0 = transformations[( byte ) *it++ ] ) == 0xFF ) {
            val = 0xFFFF;
            it++;
        } else if(( v1 = transformations[( byte ) *it++ ] ) == 0xFF ) {
            val = 0xFFFF;
        } else {
            val = ( v0 << 8 ) | v1;
        }

//        *tmp_data++ = val;
        SetUshortAtDataBlockPtr( tmp_data, val );
        ++tmp_data;
        cIdx++;
    } while( it++ != it_end );
#else
#error Incomplete implementation of adding genotype by row
#endif

}

void BasicGenotypeTable::selectMarker( uint rIdx ) {
    if( rIdx != current_dist_rIdx ) {
        resetSelectedMarker();

        current_dist_rIdx = rIdx;

        DataBlock *tmp_data = data + current_dist_rIdx * blocks_per_row;
        DataBlock val;

        ushort geno_code = 0, idx = 0;

        memset( tmp_lookup, 0, 0x10000 * sizeof( ushort ) );

        for( uint i = 0; i < blocks_per_row; ++i, ++tmp_data ) {
            if( GetUshortAtDataBlock( (val = GetDataBlockAtPtr( tmp_data )) ) == 0xFFFF ) {
                ++gt_dist.xx;
            } else {
                idx = tmp_lookup[ GetUshortAtDataBlock( val ) ];
                if( idx == 0 ) {
                    assert( geno_code != 7 );
                    if( geno_code == 3 ) {
                        // first homozygous and heterozygous found
                        assert( isGenotypeHomozygous( val ) );
                        idx = 3;
                        geno_code = 7;
                    } else if( geno_code == 5 ) {
                        // both homozygous genotypes already located
                        assert( !isGenotypeHomozygous( val ) );
                        idx = 2;
                        geno_code = 7;
                    } else if( geno_code == 1 ) {
                        // first homozygous genotype located
                        if( isGenotypeHomozygous( val ) ) {
                            idx = 3;
                            geno_code = 5;
                        } else {
                            idx = 2;
                            geno_code = 3;
                        }
                    } else if( geno_code == 2 ) {
                        assert( isGenotypeHomozygous( val ) );
                        idx = 1;
                        geno_code = 3;
                    } else if( geno_code == 0 ) {
                        if( isGenotypeHomozygous( val ) ) {
                            idx = 1;
                            geno_code = 1;
                        } else {
                            idx = 2;
                            geno_code = 2;
                        }
                    } else {
                        assert( false );
                    }
                    gt_header.header[idx] = GetUshortAtDataBlock(val);
                    tmp_lookup[ GetUshortAtDataBlock( val )] = idx;
                }
                ++gt_dist.freq[idx];
            }
        }

        row_selected = true;
    }
}

void BasicGenotypeTable::selectMarkerPair( uint rIdx1, uint rIdx2 ) {
    assert(false);
    if( !(( maIdx == rIdx2 && mbIdx == rIdx2 ) || ( maIdx == rIdx1 && mbIdx == rIdx2 ) ) ) {

    }
}

void BasicGenotypeTable::getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist ) {
    assert( false );
}

void BasicGenotypeTable::getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet &ccs, CaseControlGenotypeDistribution &ccgd ) {
    assert( false );
}

void BasicGenotypeTable::getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct ) {
    assert( false );
}

void BasicGenotypeTable::getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct ) {
    assert( false );
}

ushort BasicGenotypeTable::encodeGenotype( const string &gt ) {
    assert( gt.length() == MAX_ALLELE_COUNT );

    ushort val;
    val = ( byte ) gt[0];
    val <<= 8;
    val |= ( byte ) gt[1];

    return val;
}

const char *BasicGenotypeTable::decodeGenotype( ushort encoded_gt ) {
    if( encoded_gt != 0xFFFF && encoded_gt != 0x0000 ) {
        gt_lookup[0] = ( char )( encoded_gt >> 8 );
        gt_lookup[1] = ( char )( encoded_gt & 0xFF );
    } else {
        gt_lookup[0] = '0';
        gt_lookup[1] = '0';
    }
    return gt_lookup;
}

BasicGenotypeTable::~BasicGenotypeTable() {
    //dtor
    delete [] gt_lookup;
//    delete [] transformations;
}

}
}
