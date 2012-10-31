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
#ifndef GENOTYPERECORD_H
#define GENOTYPERECORD_H

#include <iostream>
#include <map>
#include <vector>

#include "genetics/genotype/genotype.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

class DataIterator : public iterator< input_iterator_tag, EncodedID > {
    public:
        DataIterator( EncodedID *d ) : data( d ) {}
        DataIterator( const DataIterator &i ) : data( i.data ) {}

        DataIterator &operator++() { ++data; return *this; }

        bool operator==( const DataIterator &rhs ) { return data == rhs.data; }
        bool operator!=( const DataIterator &rhs ) { return !operator==( rhs ); }

        EncodedID operator*() { return *data; }

    private:
        EncodedID *data;
};

class GenotypeDataIterator : public iterator< input_iterator_tag, EncodedID > {
    public:
        GenotypeDataIterator( GenotypeID *h, EncodedID *d ) : header( h ), cur_id( d ) {}
        GenotypeDataIterator( const GenotypeDataIterator &it ) : header( it.header ), cur_id( it.cur_id ) {}

        GenotypeDataIterator &operator++() { ++cur_id; return *this; }
        //GenotypeDataIterator &operator++(int) { GenotypeDataIterator tmp( *this ); operator++(); return tmp; }

        bool operator==( const GenotypeDataIterator &rhs ) { return cur_id == rhs.cur_id; }
        bool operator!=( const GenotypeDataIterator &rhs ) { return !operator==( rhs ); }

        GenotypeID operator*() { return header[ *cur_id ]; }

    private:
        GenotypeID *header;
        EncodedID *cur_id;
};

class GenotypeIterator : public iterator< input_iterator_tag, GenotypeID > {
    public:
        GenotypeIterator( GenotypeID *h ) : header( h ) {}
        GenotypeIterator( const GenotypeIterator &it ) : header( it.header ) {}

        GenotypeIterator &operator++() { ++header; return *this; }

        bool operator==( const GenotypeIterator &rhs ) { return header == rhs.header; }
        bool operator!=( const GenotypeIterator &rhs ) { return !operator==( rhs ); }

        GenotypeID operator*() { return *header; }

    private:
        GenotypeID *header;
};

class GenotypeRecord {
    public:
        friend class GenotypeDataIterator;

        GenotypeRecord( map< GenotypeID, EncodedID > &header, vector< EncodedID > &data );

        int getGenotypeCount() const { return data - header; }
        int getColumnCount() const { return data_end - data; }

        GenotypeIterator genotype_begin() { return GenotypeIterator( header ); }
        GenotypeIterator genotype_end() { return GenotypeIterator( data ); }

        GenotypeDataIterator genotype_data_begin() { return GenotypeDataIterator( header, data ); }
        GenotypeDataIterator genotype_data_end() { return GenotypeDataIterator( header, data_end ); }

        DataIterator begin() { return DataIterator( data ); }
        DataIterator end() { return DataIterator( data_end ); }

        virtual ~GenotypeRecord();
    protected:
        virtual void initialize( map< GenotypeID, EncodedID > &header, vector< EncodedID > &data );

        char *mem_block;
        GenotypeID *header;
        EncodedID *data, *data_end;
    private:
};

}
}

#endif // GENOTYPERECORD_H
