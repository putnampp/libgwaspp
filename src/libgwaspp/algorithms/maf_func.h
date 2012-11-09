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
#ifndef MAF_FUNC_H
#define MAF_FUNC_H

#include <fstream>

#include "genetics/genetic_data.h"
#include "genetics/genotype/geno_table.h"
#include "util/time/timing.h"

#include "algorithms/computation_engine.h"

namespace libgwaspp {
namespace algorithms {

using namespace libgwaspp::genetics;
using namespace std;

inline void maf( const frequency_table & ft, double & tot, double & maf ) {
    tot = ft.aa;
    maf = 2.0 * tot;
    tot += ft.ab;
    maf += ft.ab;
    tot += ft.bb;
    maf /= tot;
    if( maf < 0.5 ) { maf = 1.0 - maf; }
}

const int POSSIBLE_ENC = 1 << ( 8 * sizeof( ushort ) );

void computeMAF( GenoTable &gt );
void computeMAF( GeneticData *gd, const set<string> & row_ids, const set<string> & column_ids );

void maf_all( void *input, void *output );
void maf_from_distribution( void *input, void *output );

void maf( void *input, void *output );
void maf_all( IndexedInput &input, void *output );

void compute_maf_perform( GeneticData *gd, ostream *out);

void select_cc_maf( GeneticData *gd, ostream *out);
void inline_cc_maf( GeneticData *gd, ostream *out);

void inline_maf_print( GeneticData *gd, ostream *out );

}
}

#endif // MAF_FUNC_H
