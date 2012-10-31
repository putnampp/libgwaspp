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
#ifndef SERIALIZABLE_H
#define SERIALIZABLE_H

#include <fstream>
#include <iostream>

namespace libgwaspp {


/**
    Interface which can be utilized for making objects serializable

    Every child class must have a save and load function
    @author Partick Putnam
*/
struct loadable {
    virtual std::istream &load( std::istream &is ) = 0;
};

struct saveable {
    virtual std::ostream &save( std::ostream &os ) = 0;
};

struct builder {
    virtual std::istream &build( std::istream &is ) = 0;
};

template< typename O >
class streambuilder {
    public:
        virtual std::istream &build( std::istream &is ) = 0;
        virtual O *getObject() = 0;
    protected:
        O * obj;
};

template< typename C >
struct col_data : public streambuilder< C > {
    virtual std::istream& parseHeader( std::istream& is ) = 0;
};

template< typename R, typename C, typename D >
class TabularRecordParser {
    public:
        TabularRecordParser( col_data< R > * r, col_data< C > * c, streambuilder< D > * d ) : rows( r ), columns( c ), data( d ) {}

        bool parseRecords( std::istream& is);

        virtual ~TabularRecordParser() {}
    protected:
        virtual std::istream &parseColumnHeader( std::istream &is ) = 0;
        virtual std::istream &parseRowHeader( std::istream &is ) = 0;
        virtual std::istream &parseRowData( std::istream &is ) = 0;
        virtual std::istream &parseData( std::istream &is ) = 0;
        streambuilder< R > * rows;
        streambuilder< C > * columns;
        streambuilder< D > * data;

        istringstream iss;
        string line;
};

bool TabularRecordParser::parseRecords( std::istream &is ) {
    assert( rows != NULL && columns != NULL && data != NULL);

    if( !is.good() ) {
        cout << "Unable to open Phenotype File" << endl;
        return false;
    }

    getline( is, line);

    iss.str( line );
    iss.clear();

    parseRowHeader( iss );
    while( !iss.eof() ) {
        parseColumnHeader( iss );
    }

    while( !ifs.eof() ) {
        getline(ifs, line);
        iss.str( line );
        iss.clear();
        parseRowData( iss );

        while( !iss.eof() ) {
            parseData(iss);
        }
    }

    ifs.close();
    return true;
}

struct serializable : public loadable, public saveable { };

}

#endif // SERIALIZABLE_H
