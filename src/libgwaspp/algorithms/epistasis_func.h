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
#ifndef EPISTASIS_FUNC_H
#define EPISTASIS_FUNC_H

#include <cstring>
#include <cmath>
#include <fstream>

#include "genetics/genetic_data.h"
#include "genetics/genotype/common_genotype.h"
#include "genetics/genotype/common_genotype_func.h"
#include "genetics/analyzable/case_control_set.h"

#include "util/time/timing.h"

#include "algorithms/computation_engine.h"

#include "boost/format.hpp"

//#define MATHLIB_STANDALONE
//#include "Rmath.h"

namespace libgwaspp {
namespace algorithms {

using namespace libgwaspp::genetics;

//struct marginal_information {
//    frequency_table margins, cases, controls;
//    double dMarginalEntropy, dMarginalEntropy_Y;
//};
typedef pair< uint, uint > SNPPair;
typedef pair< SNPPair, double> SNPInteractionPair;

const int GT_COUNT = 3;
const int GT_BUFFER_COUNT = 2;
const int GT_BUFFER_SIZE = GT_COUNT *GT_BUFFER_COUNT;

void epistasis_all( void *input, void *output );
void epistasis_all2( void *input, void *output );

void ContingencyDebug( void * input, void * output );
void ContingencyPerformance( void * input, void * output );

void EpistasisDebug( void * input, void * output );
void EpistasisPerformance( void * input, void * output );

void computeGTest( GenoTable & gt, marginal_information * pMargins, uint nIndivids, vector< SNPInteractionPair > & passingThreshold, vector< double > & zval );

void computeMargins( GenoTable & gt, int nIndivids, marginal_information *& pMargins, int & nMarkerCount );
void computeBoost( GeneticData *gd, ostream *out );

double pairwise_epi_test( const CONTIN_TABLE_T &_case, const CONTIN_TABLE_T &_ctrl);
double pairwise_epi_test( const CONTIN_TABLE_T &_case, const CONTIN_TABLE_T &_ctrl, const marginal_information & m1, const marginal_information & m2 );

}
}

#endif // EPISTASIS_FUNC_H
