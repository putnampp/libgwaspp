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
#ifndef COMPUTATIONENGINE_H
#define COMPUTATIONENGINE_H

#include <iostream>
#include <fstream>
#include <memory>
#include <set>

#include "genetics/analyzable/case_control_set.h"
#include "genetics/genetic_data.h"
#include "util/time/timing.h"

using namespace std;
using namespace libgwaspp::genetics;

namespace libgwaspp {
namespace algorithms {

struct BasicInput {
    GeneticData *gd;
    set<string> * marker_ids, * individual_ids;

    BasicInput() : gd( NULL ), marker_ids( new set<string>() ), individual_ids( new set<string>() ) {}
    BasicInput( GeneticData *g, set<string> * m_ids, set<string> * i_ids ) : gd( g ), marker_ids( m_ids ), individual_ids( i_ids ) {}
    BasicInput( BasicInput *bi ) : gd( bi->gd ), marker_ids( bi->marker_ids ), individual_ids( bi->individual_ids ) {}

    virtual ~BasicInput() { marker_ids->clear(); individual_ids->clear(); }
};

struct IndexedInput : BasicInput {
    auto_ptr< set<int> > midx, iidx;

    IndexedInput( BasicInput *bi ) : BasicInput( bi ), midx( new set<int>() ), iidx( new set<int>() ) {}

    ~IndexedInput() { midx->clear(); iidx->clear(); }
};

struct CaseControlInput : BasicInput {
    CaseControlSet ccs;

    CaseControlInput( BasicInput *bi ) : BasicInput( bi ), ccs( gd->getGenotypeTable()->getColumnSet() ) {}
};

void compute( void ( *f )( void *, void * ), void *input, void *output );
void compute( void ( *f )( GeneticData *, ostream *), GeneticData *gd, ostream * out );

}
}

#endif // COMPUTATIONENGINE_H
