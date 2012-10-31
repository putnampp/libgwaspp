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
#ifndef GENOTYPETABLE2_H
#define GENOTYPETABLE2_H

#include <iostream>
#include <iterator>
#include <cstring>
#include <cmath>
#include <cassert>

#include "libgwaspp.h"
#include "genetics/genotype/genotype.h"
#include "util/table/table.h"
#include "util/index_set/indexer.h"

using namespace util;
using namespace std;

namespace libgwaspp {
namespace genetics {


class GenoIterator : public iterator< input_iterator_tag, byte > {
    public:
        GenoIterator( byte *p, int step = 1 ) : d( p ), _step( step ) {}
        GenoIterator( const GenoIterator &it ) : d( it.d ) {}

        GenoIterator &operator++() { d += _step; return *this; }
        //GenoIterator &operator++( int ) { GenoIterator tmp( *this ); operator++(); return tmp; }

        bool operator==( const GenoIterator &rhs ) { return d == rhs.d; }
        bool operator!=( const GenoIterator &rhs ) { return !operator==( rhs ); }

        byte operator*() { return *d; }

    private:
        byte *d;
        int _step;
};


class GenotypeTable2 : public Table< byte > {
    public:
        GenotypeTable2( indexer *markers, indexer *individuals ) : Table< byte >( markers, individuals ), headers( NULL ), gt_lookup(NULL) {
            initialize();
            cout << "Allocated Genotype Table" << endl;
        }

        byte operator()( int rIdx, int cIdx );
//        byte operator()( const string &row_id, const string &column_id );

        void addGenotype( int rIdx, int cIdx, const string &gt );

        void addGenotypeRow( int rIdx, string::const_iterator & it, string::const_iterator & it_end, char delim);

        byte encodeGenotype( const string &gt );
        const char *decodeGenotype( byte encoded_gt );

//        GenoIterator row_start( int rIdx ) { return GenoIterator( data + rIdx * max_column ); }
//        GenoIterator row_end( int rIdx ) { return GenoIterator( data + ( rIdx + 1 ) * max_column ); }
//
//        GenoIterator row_start( const string &rID ) { return row_start( (*rows)( rID ) ); }
//        GenoIterator row_end( const string &rID ) { return row_end( (*rows)( rID ) ); }
//
//        GenoIterator column_start( int cIdx ) { return GenoIterator( data + cIdx, max_column ); }
//        GenoIterator column_end( int cIdx ) { return GenoIterator( data + ( max_row - 1 ) * max_column + cIdx, max_column ); }
//
//        GenoIterator column_start( const string &cID ) { return column_start( (*columns)( cID ) ); }
//        GenoIterator column_end( const string &cID ) { return column_end( (*columns)( cID ) ) ; }



        virtual ~GenotypeTable2();
    protected:
        void initialize();
        void setupAlphabetTranslation();

        byte *headers;
        byte alphabet_translation[256];
//        char **gt_lookup;
        char * gt_lookup;

        static const string alphabet;
        static const string err_alphabet;

        int alphabet_size, bits_per_gt, bytes_per_gt;
        int header_bit_len, header_byte_len;
        int bits_per_row, bytes_per_row;
        int bits_per_data;
        int possible_genotypes_size, local_genotype_count;

        int data_size, header_size, lookup_size, gt_size;

    private:
        //GenotypeTable2( int r, int c) {}
};

}
}

#endif // GENOTYPETABLE2_H
