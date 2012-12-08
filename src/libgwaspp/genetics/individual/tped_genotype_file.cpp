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
#include "genetics/individual/tped_genotype_file.h"

namespace libgwaspp {
namespace genetics {

TpedGenotypeFile::TpedGenotypeFile() : buffer(NULL)
{
    //ctor
}

void TpedGenotypeFile::guessExpectedSizes( istream *iFile, int &mcnt, int &gt_width, char delim ) {
    iFile->seekg( 0, ios::beg );

    getline( *iFile, line );
    long len = iFile->tellg();
    getline( *iFile, line );
    iFile->seekg( 0, ios::end );
    len = iFile->tellg() - len;

    parser.str( line );
    parser.clear();

    for( int i = 0; i < 4; i++ ) { getline( parser, tok, delim ); }

    mcnt = ( int )( len / line.length() ) + 30;

    getline( parser, tok, delim );
    gt_width = tok.length();
}

bool TpedGenotypeFile::parseHeader( istream *iFile, GeneticData *gd, char delim ) {
//    if( delim == 0 ) return false;
//    bool success = true;
//
//    iFile->seekg( 0, ios::beg );
//    getline( *iFile, line );
//
//    // skip the first 4 columns
//    // assumed file format:
//    // chrom\trs\tgenetic distance\tbase-pair position
//
//    string::const_iterator s_it = line.begin(), s_end = line.end();
//    int delim_count = 0;
//
//    while(delim_count < 4 ) {
//        switch( *s_it++) {
//            case ' ':
//            case '\t':
//                ++delim_count;
//                break;
//            default:
//                break;
//        }
//    }
//
//
//    individ_indexes.clear();
//    // remaining columns correspond to individual ids
//    int idx = 0;
//    while( s_it < s_end ) {
//        s_it += 4;  // tped alleles are assumed to be bi-allelic and separated
//
//        individ_indexes.push_back( idx++ );
//    }
//
//    gd->createGenotypedIndividuals( individ_indexes );
//    individ_indexes.clear();
//
//    iFile->seekg( 0, ios::beg );
//    return success;
    return true;
}

bool TpedGenotypeFile::parseNextMarkerRecord( istream *iFile, GeneticData *gd, char delim ) {
    if( delim == 0 ) return false;
    bool success = true;

    getline( *iFile, line );
    parser.str( line );
    parser.clear();

    // Parse each marker information and build the marker
    string rs, chr, strand, allele = "";
    uint pos;

    getline( parser, chr, delim );
    getline( parser, rs, delim );
    getline( parser, strand, delim );
    parser >> pos;
    getline( parser, tok, delim );  // read delim following position

    int row_idx = gd->getMarkerIndex( gd->createMarker( rs, chr, ( uint )pos, pos + 1, -1.0, allele ) );

    marker_indexes.push_back( row_idx );

    return success;
}

bool TpedGenotypeFile::parseNextGenotypeRecord( istream *iFile, GeneticData *gd, char delim ) {
    if( delim == 0 ) return false;
    bool success = true;

    getline( *iFile, line );

    trim(line);

    string::const_iterator it = line.begin();
    string::const_iterator it_end = line.end();
    int delim_count = 0;
    while( delim_count < 4 ) {
        //if( *it++ == delim ) { ++delim_count; }
        switch( *it++ ) {
            case ' ':
            case '\t':
                ++delim_count;
                break;
            default:
                break;
        }
    }

    if( buffer == NULL ) {
        buffer_size = (((it_end - it) + 1) >> 2) * 3;
        buffer = new char[ buffer_size + 1 ];
        memset( buffer, 0, buffer_size + 1 );
    }

    tmp_buffer = buffer;

    bool set_delim = false;
    // remove separating delimiter between alleles
    while( it < it_end ) {
        switch( *it ) {
            case '1':
                *tmp_buffer++ = 'A';
                break;
            case '2':
                *tmp_buffer++ = 'C';
                break;
            case '3':
                *tmp_buffer++ = 'G';
                break;
            case '4':
                *tmp_buffer++ = 'T';
                break;
            default:
                *tmp_buffer++ = *it;
                break;
        }
        it += 2;
        if( set_delim ) {
            *tmp_buffer++ = delim;
            set_delim = false;
        } else {
            set_delim = true;
        }
    }

    gd->addGenotypeRow( row_idx++, buffer, buffer + buffer_size, delim );

    return success;
}

TpedGenotypeFile::~TpedGenotypeFile()
{
    //dtor
    delete [] buffer;
}

}
}
