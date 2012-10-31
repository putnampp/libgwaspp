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
#include "genetics/individual/individual_phenotype_file.h"

namespace libgwaspp {
namespace genetics {

bool IndividualPhenotypeFile::populateGeneticData( string &filename, GeneticData *gd, char delim ) {
    bool success = true, is_gz = false;
    istream *iFile;
    string ext = filename.substr(( int )filename.find_last_of( "." ) + 1 );
    if( ext.compare( "gz" ) == 0 || ext.compare( "GZ" ) == 0 ) {
        iFile = new igzstream( filename.c_str() );
        is_gz = true;
    } else {
        iFile = new ifstream( filename.c_str() );
    }

    if( parseHeader( iFile, gd, delim ) ) {
        while( !iFile->eof() ) {
            if( ! parseNextRecord( iFile, gd, delim ) ) {
                success = false;
                break;
            }
            iFile->peek();
        }
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

bool IndividualPhenotypeFile::parseHeader( istream *iFile, GeneticData *gd, char delim ) {
    iFile->seekg( 0, ios::beg );
    getline( *iFile, line );

    parser.str( line );
    parser.clear();

    // skip the first 6 columns of the header
    // assumed to be:
    // "idx"\t"id"\t"sex"\t"age"\t"fid"\t"mid"
    for( int i = 0; i < 6; i++ ) { getline( parser, tok, delim ); }

    // remaining variable number of columns represent traits which have been phenotyped
    while( getline( parser, tok, delim ) ) {
        cout << "Header token: " << tok << endl;

        gd->addOrGetPhenotype( tok );
    }

    return true;
}

bool IndividualPhenotypeFile::parseNextRecord( istream *iFile, GeneticData *gd, char delim ) {
    bool success = true;
    //string line;

    getline( *iFile, line );

    //istringstream parser( line );
    parser.str( line );
    parser.clear();

//    PhenotypeTree *ptree = individs->getPhenotypeTree();

    //float age;

    getline( parser, tok, delim );  // skip index column
    getline( parser, id, delim );   // get the current individuals id
    getline( parser, sex, delim );  // get the current individuals sex
    //getline( parser, age, '\t');
    //parser >> std::dec >> age; // get the current individuals age
    getline( parser, age, delim );  // skip the \t after age
    getline( parser, fid, delim );  // get the current individuals fid
    getline( parser, mid, delim );  // get the current individuals mid

    Individual *ind, *mInd, *fInd;
    ind = gd->findOrCreateIndividual( id );
    mInd = gd->findOrCreateIndividual( mid );
    fInd = gd->findOrCreateIndividual( fid );

//    cout << id << " <- " << fid << " + " << mid << endl;

    ind->setMother( mInd );
    ind->setFather( fInd );

    ind->setAge( age );
    ind->setSex( sex );

    vector<PhenotypeNode *>::iterator n_it = gd->phenotype_begin();
    const PhenotypeValue *val;

    while( getline( parser, tok, delim ) ) {
        val = gd->addOrGetPhenotypeValue( *n_it, tok );
        ind->addPhenotype( val );
        n_it++;
    }

    return success;
}

IndividualPhenotypeFile::~IndividualPhenotypeFile() {
    line.clear();
    tok.clear();
    id.clear();
    sex.clear();
    mid.clear();
    fid.clear();
}

}
}
