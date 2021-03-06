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
#ifndef GENO_TABLE_H
#define GENO_TABLE_H


#include <iostream>

#include "libgwaspp.h"
#include "util/table/table.h"
#include "util/index_set/indexer.h"

#include "genetics/genotype/genotype.h"
#include "genetics/genotype/pairwise_marker_analyzable.h"
#include "genetics/genotype/single_marker_analyzable.h"
#include "genetics/genotype/case_control_selectable.h"

namespace libgwaspp {
namespace genetics {

class GenoTable : public util::Table < DataBlock >, public PairwiseMarkerAnalyzable, public SingleMarkerAnalyzable, public virtual CaseControlSelectable {
public:
    GenoTable( util::indexer *markers, util::indexer *individuals ) : util::Table< DataBlock >( markers, individuals ), data_per_block( 0 ), blocks_per_row( 0 ), total_block_count( 0 ), bytes_per_row( 0 ) {}

    virtual void addGenotype( int rIdx, int cIdx, const string &gt ) = 0;
    virtual void addGenotypeRow( int rIdx, string::const_iterator &it, string::const_iterator &it_end, char delim ) = 0;
    virtual void addGenotypeRow( int rIdx, const char *p_begin, const char *p_end, char delim ) = 0;

    virtual ushort encodeGenotype( const string &gt ) = 0;
    virtual const char *decodeGenotype( ushort encoded_gt ) = 0;

    inline const char *getCallAt( uint marker_idx, uint individ_idx ) { return decodeGenotype( getGenotypeAt( marker_idx, individ_idx )); }
    DataBlock getGenotypeAt( uint marker_idx, uint individ_idx ) { return (*this)( marker_idx, individ_idx ); }

    uint getPossibleGenotypeCount() const { return possible_genotypes_size; }
    virtual ushort *possible_genotypes_begin() const { return beg; }
    virtual ushort *possible_genotypes_end() const { return end; }

    virtual bool isGenotypeHomozygous( ushort enc ) = 0;

    virtual ~GenoTable() {
        delete [] beg;
    }

protected:
    byte transformations[256];

    ulong data_per_block, blocks_per_row, total_block_count;
    uint bytes_per_row;
    uint bits_per_data;
    uint alphabet_size, possible_genotypes_size;
    ulong data_size;

    ushort *beg, *end;

    static const int BITS_PER_BLOCK;
    static const int BYTES_PER_BLOCK;

    static const std::string alphabet;
    static const std::string err_alphabet;
};

}
}

#endif // GENO_TABLE_H
