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
#include "genetics/individual/individual_genotype_file.h"

namespace libgwaspp {
namespace genetics {

IndividualGenotypeFile::IndividualGenotypeFile() : GeneticDataFile(),  row_idx( 0 ), col_idx( 0 ) {
    //ctor
}

bool IndividualGenotypeFile::populateGeneticData( string &filename, GeneticData *gd, char delim ) {
    bool success = true, is_gz = false;
    istream *iFile;
    string ext = filename.substr(( int )filename.find_last_of( "." ) + 1 );
    if( ext.compare( "gz" ) == 0 || ext.compare( "GZ" ) == 0 ) {
        iFile = new igzstream( filename.c_str() );
        is_gz = true;
    } else {
        iFile = new ifstream( filename.c_str() );
    }

    INIT_LAPSE_TIME;
    RECORD_START;

    if( parseHeader( iFile, gd, delim ) ) {

        long first_rec = iFile->tellg();

        RECORD_STOP;
        PRINT_LAPSE(cout, "Parse Individuals Time Lapse: " );
        cout << endl;
        // build markers
        marker_indexes.clear();

        RECORD_START;
        while( !iFile->eof() ) {
            if( ! parseNextMarkerRecord( iFile, gd, delim ) ) {
                success = false;
                break;
            }
            iFile->peek();
        }

        assert(( int )marker_indexes.size() == gd->getMarkerCount() );

        RECORD_STOP;
        PRINT_LAPSE(cout, "Parse Marker Time Lapse: " );
        cout << endl;

        RECORD_START;
        gd->createGenotypedMarkers( marker_indexes );

        assert( gd->getGenotypedMarkersCount() == gd->getMarkerCount() );

        RECORD_STOP;
        PRINT_LAPSE(cout, "Create Genotyped Marker Time Lapse: " );
        cout << endl;

        marker_indexes.clear();

        RECORD_START;
        gd->updateGenotypeTable();
        RECORD_STOP;
        PRINT_LAPSE(cout, "Lapsed time for updating Genotype Table for genotyped individuals and markers: " );
        cout << endl;

        iFile->seekg( first_rec, ios::beg );
        iFile->clear();

        RECORD_START;
        // build genotypes
        row_idx = 0;
        while( !iFile->eof() ) {
            if( ! parseNextGenotypeRecord( iFile, gd, delim ) ) {
                success = false;
                break;
            }
            iFile->peek();
        }

        RECORD_STOP;
        PRINT_LAPSE(cout, "Lapsed time for compacting Genotype Data: " );
        cout << endl;
    } else {
        success = false;
    }

    if( is_gz ) {
        (( igzstream * )iFile )->close();
    } else {
        (( ifstream * )iFile )->close();
    }

    delete iFile;

    return success;
}

void IndividualGenotypeFile::guessExpectedSizes( istream *iFile, int &mcnt, int &gt_width, char delim ) {
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

bool IndividualGenotypeFile::parseHeader( istream *iFile, GeneticData *gd, char delim ) {
    if( delim == 0 ) return false;
    bool success = true;

    iFile->seekg( 0, ios::beg );
    getline( *iFile, line );

    parser.str( line );
    parser.clear();

    // skip the first 4 columns
    // assumed file format:
    // rs\tchrom\tpos\tstrand
    for( int i = 0; i < 4; i++ ) { getline( parser, tok, delim ); }

    individ_indexes.clear();
    // remaining columns correspond to individual ids
    int idx = 0;
    while( getline( parser, tok, delim ) ) {
        idx = gd->findOrCreateIndividualIndex( tok );
        individ_indexes.push_back( idx );
    }

    gd->createGenotypedIndividuals( individ_indexes );
    individ_indexes.clear();
    return success;
}

bool IndividualGenotypeFile::parseNextMarkerRecord( istream *iFile, GeneticData *gd, char delim ) {
    if( delim == 0 ) return false;

    bool success = true;

    getline( *iFile, line );

    parser.str( line );
    parser.clear();

    // Parse each marker information and build the marker
    string rs, chr, strand, allele = "";
    uint pos;

    getline( parser, rs, delim );
    getline( parser, chr, delim );
    parser >> pos;
    getline( parser, tok, delim );  // read tab following position
    getline( parser, strand, delim );

    int row_idx = gd->getMarkerIndex( gd->createMarker( rs, chr, ( uint )pos, pos + 1, -1.0, allele ) );

    marker_indexes.push_back( row_idx );

    return success;
}

bool IndividualGenotypeFile::parseNextGenotypeRecord( istream *iFile, GeneticData *gd, char delim ) {
    if( delim == 0 ) return false;

    bool success = true;

    getline( *iFile, line );

    trim(line);

    string::const_iterator it = line.begin();
    string::const_iterator it_end = line.end();
    int data_start = 0, delim_count = 0;
    while( delim_count < 4 ) {
        if( *it == delim ) { ++delim_count; }
        ++it;
        ++data_start;
    }

    gd->addGenotypeRow( row_idx++, it, it_end, delim );

    return success;
}

IndividualGenotypeFile::~IndividualGenotypeFile() {
    //dtor
    marker_indexes.clear();
    individ_indexes.clear();
}

}
}
