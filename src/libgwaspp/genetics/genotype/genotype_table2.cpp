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
#include "genetics/genotype/genotype_table2.h"

namespace libgwaspp {
namespace genetics {

const string GenotypeTable2::alphabet = GENOTYPE_ALPHABET;
const string GenotypeTable2::err_alphabet = "0?N";

byte GenotypeTable2::operator()( int rIdx, int cIdx ) {
    byte v = get_gt_value_at( data + rIdx * bytes_per_row, cIdx );

    if( v == 0 ) return 0;
    --v;
    return *( headers + rIdx * header_byte_len + v );
}

void GenotypeTable2::addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim ) {
#if MAX_ALLELE_COUNT == 2
    byte enc = 0, tmp_e = 0;
    char c1, c2;
    byte val1 = 0, val2 = 0;

    int max_set_count = 1;
    char enc_set[ 256 ];
    memset( enc_set, ( char ) - 1, 256 );

    byte *tmp_data = data + rIdx * bytes_per_row;

    // make sure the data row is "unknown"
    memset( tmp_data, 0, bytes_per_row );

//    int cIdx = 0;
//    do {
//        c1 = *it++;
//        c2 = *it++;
//
//        // no point in respecifying that an "unknown" column is "unknown", therefore only modify columns which are "known"
//        if( !( ( val1 = alphabet_translation[( byte )c1] ) == 255 || ( val2 = alphabet_translation[( byte )c2] ) == 255) ) {
//            enc = val1 * alphabet_size + val2 + 1;
//            tmp_e = enc_set[enc];
//            if(( tmp_e = enc_set[ enc ] ) == 255 ) {
//                assert( max_set_count < 4 );
//                tmp_e = max_set_count++;
//                enc_set[ enc ] = tmp_e;
//            }
//
//            assert( tmp_e != 0 );
//            set_gt_value_at( tmp_data, cIdx, tmp_e );
//        }
//        cIdx++;
//    } while( it++ != it_end );


    // method offers fewer pointer dereferences
    byte val = 0, offset = 0;

    do {
        c1 = *it++;
        c2 = *it++;

        if( offset >= 8 ) {
            *tmp_data = val;
            ++tmp_data;
            offset = 0;
            val = 0;
        }

        // no point in respecifying that an "unknown" column is "unknown", therefore only modify columns which are "known"
        if( !( ( val1 = alphabet_translation[( byte )c1] ) == 255 || ( val2 = alphabet_translation[( byte )c2] ) == 255) ) {
            enc = val1 * alphabet_size + val2 + 1;
            tmp_e = enc_set[enc];
            if(( tmp_e = enc_set[ enc ] ) == 255 ) {
                assert( max_set_count < 4 );
                tmp_e = max_set_count++;
                enc_set[ enc ] = tmp_e;
            }

            assert( tmp_e != 0 );
            val |= (tmp_e << offset);
        }

        offset += 2;
    } while( it++ != it_end );

    if( offset > 0 ) {
        *tmp_data = val;
    }

    // construct header
    byte *tmp_header = headers + rIdx * header_byte_len;
    for( int i = 0, j = 1; i < 256 && j < max_set_count; ++i ) {
        if(( tmp_e = enc_set[ i ] ) != 255 ) {
            --tmp_e;
            *( tmp_header + tmp_e ) = ( byte ) i;
        }
    }
#else
#error Incomplete implementation of adding genotype by row
#endif
}

void GenotypeTable2::addGenotype( int rIdx, int cIdx, const string &gt ) {
    byte tmp_e = 0;
    byte enc = encodeGenotype( gt );

    if( enc != 0 ) {
        byte *tmp_head = headers + rIdx * header_byte_len;
        bool found = false;
        for( int i = 1; i < header_byte_len && !found; ++i, ++tmp_head ) {
            if( *tmp_head == enc ) {
                tmp_e = i;
                found = true;
            } else if( *tmp_head == 0 ) {
                tmp_e = i;
                *tmp_head = enc;
                found = true;
            }
        }

        assert( found );
    }

    set_gt_value_at( data + rIdx * bytes_per_row, cIdx, tmp_e );
}

byte GenotypeTable2::encodeGenotype( const string &gt ) {
#if MAX_ALLELE_COUNT == 2
    assert( gt.length() == MAX_ALLELE_COUNT );
    int idx0 = alphabet_translation[( byte ) gt[0] ],
               idx1 = alphabet_translation[( byte ) gt[1] ];

    if( idx0 == -1 || idx1 == -1 ) return 0;

    return idx0 * alphabet_size + idx1 + 1;
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

const char *GenotypeTable2::decodeGenotype( byte encoded_gt ) {
//    return gt_lookup[encoded_gt];
    return ( gt_lookup + encoded_gt * gt_size );
}

void GenotypeTable2::initialize() {
    cout << "Initializing GenotypeTable2 ... " << endl;
    alphabet_size = ( int ) alphabet.length();
    possible_genotypes_size = pow( alphabet_size, ( double ) MAX_ALLELE_COUNT ) + 1;   // +1 for 0 -> "unknown" genotype

    cout << "Alphabet Size: " << alphabet_size << endl;
    cout << "Possible Genotype Size: " << possible_genotypes_size << endl;

    assert( possible_genotypes_size < 256 );

    bits_per_gt = ( int ) ceil( log(( double )possible_genotypes_size ) / log( 2.0 ) );

    cout << "Bits per Genotype: " << bits_per_gt << endl;

    local_genotype_count = ( int ) GENOTYPE_FORM_COUNT - 1; // -1 for common 0 across all local
    bits_per_data = ( int ) ceil( log(( double ) local_genotype_count ) / log( 2.0 ) );

    cout << "Bits per Data: " << bits_per_data << endl;

    bits_per_row = bits_per_data * max_column;
    bytes_per_row = byte_count( bits_per_row );

    cout << "Bits per Data Row: " << bits_per_row << endl;
    cout << "Bytes per Data Row: " << bytes_per_row << endl;

    header_bit_len = local_genotype_count * bits_per_gt;  // Each Genotype Record can have at most 3 unique encodings (AA, AB, BB) + 1 unknown (00)
    header_byte_len = local_genotype_count;

    if( data != NULL ) {
        delete [] data;
    }
    data_size = max_row * bytes_per_row;

    cout << "Total Table size: " << data_size << " (bytes)" << endl;

    data = new byte[ data_size ];
    memset( data, 0, data_size );

    if( headers != NULL ) {
        delete [] headers;
    }
    header_size = max_row * header_byte_len;

    cout << "Total Header size: " << header_size << " (bytes)" << endl;

    headers = new byte[ header_size ];
    memset( headers, 0, header_size );

    if( gt_lookup != NULL ) {
        delete [] gt_lookup;
    }

    gt_size =  MAX_ALLELE_COUNT + 1;
    int idx, offset;
    lookup_size = gt_size * possible_genotypes_size;

    gt_lookup = new char[ lookup_size ];    // allocate space for all possible genotypes as null-terminated character sequences
    memset( gt_lookup, ( char )0, lookup_size );
    for( int i = 0; i < possible_genotypes_size; ++i ) {
        if( i > 0 ) {
            idx = i - 1;
            for( int k = 0, l = MAX_ALLELE_COUNT - 1; k < MAX_ALLELE_COUNT; ++k, --l ) {
                offset = idx % alphabet_size;
                gt_lookup[ i *gt_size + l ] = alphabet[offset];
                idx -= offset;
                idx /= alphabet_size;
            }
        } else {
            for( int k = 0; k < MAX_ALLELE_COUNT; ++k ) {
                gt_lookup[ k ] = '0';
            }
        }

        cout << ( gt_lookup + i * gt_size ) << " -> " << i << endl;
    }
    cout << "Initialized Genotype Lookups: " << lookup_size << " (bytes)" << endl;

    setupAlphabetTranslation();
}

void GenotypeTable2::setupAlphabetTranslation() {
    memset( alphabet_translation, ( char ) - 1, 256 );

    int idx = 0;
    for( string::const_iterator it = alphabet.begin(); it != alphabet.end(); ++it, ++idx ) {
        alphabet_translation[( byte ) *it ] = ( byte )idx;
    }
}

GenotypeTable2::~GenotypeTable2() {
    delete [] headers;
    delete [] gt_lookup;
}

}
}
