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
#include "genetics/individual/tfam_annotation_file.h"

namespace libgwaspp {
namespace genetics {

TFamAnnotationFile::TFamAnnotationFile() : GeneticDataFile() {
    //ctor
}

bool TFamAnnotationFile::populateGeneticData( string &filename, GeneticData   *gd, char delim ) {
    if( filename.empty() ) {
        gd->setCaseControlSet(NULL, NULL);
        return true;
    }

    ifstream iFile( filename.c_str() );

    set<string> cases, controls;

    delim = ' ';
    string id;
    while( !iFile.eof() ) {
        getline( iFile, line );
        if( line == "" )
            continue;
        parser.str( line );
        parser.clear();

        getline( parser, id, delim );  // family id
        getline( parser, tok, delim );   // individual id
        id += "-" + tok;
        getline( parser, tok, delim);   // paternal id
        getline( parser, tok, delim);   // maternal id
        getline( parser, tok, delim);   // sex
        getline( parser, tok, delim);   // disease?

        if( id == "-1" ) {
            cout << line << endl;
        }

        switch( tok[0] ) {
        case '1':
            cases.insert( id );
            break;
        case '0':
            controls.insert( id );
            break;
        default:
            break;
        }
    }

    cout << "Setting " << cases.size() << " cases." << endl;
    cout << "Setting " << controls.size() << " controls." << endl;
    gd->setCaseControlSet( &cases, &controls);

    cases.clear();
    controls.clear();

    iFile.close();
    return true;
}

TFamAnnotationFile::~TFamAnnotationFile() {
    //dtor
}

}
}
