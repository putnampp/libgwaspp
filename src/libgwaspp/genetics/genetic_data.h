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
#ifndef GENETICDATA_H
#define GENETICDATA_H

#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cstring>
#include <iterator>
#include <cstring>
#include <cstdlib>

#include "libgwaspp.h"

#include "genetics/individual/individual_collection.h"
#include "genetics/genotype/genotype_record.h"
#include "genetics/genotype/genotype_record_factory.h"
#include "genetics/genotype/genotype_collection.h"
#include "genetics/marker/marker_collection.h"
#include "genetics/genotype/genotype_table2.h"
#include "util/index_set/indexer.h"

#include "genetics/genotype/geno_table.h"
#include "genetics/analyzable/case_control_set.h"

#if COMPRESSION_LEVEL == 1
#include "genetics/genotype/compressed_genotype_table.h"
#elif COMPRESSION_LEVEL == 2
#include "genetics/genotype/compressed_genotype_table2.h"
#elif COMPRESSION_LEVEL == 3
#include "genetics/genotype/compressed_genotype_table3.h"
#else
#include "genetics/genotype/basic_genotype_table.h"
#endif

#include "algorithms/genetic_data_func.h"

using namespace std;
using namespace util;
using namespace libgwaspp::algorithms;

namespace libgwaspp {
namespace genetics {

class GeneticData {
    public:
        GeneticData( bool isPhased = false );

        int getIndividualCount() const { return individuals->getCount(); }
        int getPhenotypeCount() const { return individuals->getPhenotypeCount(); }
        int getMarkerCount() const { return markers->getMarkerCount(); }

        int getGenotypedIndividualsCount() const { return genotyped_individs->included_size(); }
        int getGenotypedMarkersCount() const { return genotyped_markers->included_size(); }

        int getGenotypeCount() const { return getMarkerCount() * getIndividualCount(); }

        bool isPhased() const { return phased; }

        void setExpGenotypeWidth( int mwidth );
        void setExpMarkerCount( int mcnt );

        void addGenotypeRow( int r, string::const_iterator & it, string::const_iterator & it_end, char delim) { geno_tbl->addGenotypeRow( r, it, it_end, delim); }
        void addGenotypeRow( int rIdx, const char * p_begin, const char * p_end, char delim ) { geno_tbl->addGenotypeRow( rIdx, p_begin, p_end, delim); }

        void addGenotype( int r, int c, const string & s );
        const char * getGenotype( int r, int c );

        // IndividualCollection Wrapper function
        Individual *findOrCreateIndividual( string &id ) { return individuals->findOrCreateIndividual( id ); }
        Individual *findIndividualByIndex( int idx ) { return individuals->getIndividualAt( idx ); }
        int findOrCreateIndividualIndex( const string & id ) { return (*individuals)(id); }

        // Phenotype Tree Wrapper functions
        PhenotypeNode *addOrGetPhenotype( PhenotypeID &pheno_id ) { return individuals->getPhenotypeTree()->addOrGetPhenotype( pheno_id ); }
        const PhenotypeValue *addOrGetPhenotypeValue( PhenotypeNode *node, PhenotypeValue &pheno_val ) { return individuals->getPhenotypeTree()->addOrGetPhenotypeValue( node, pheno_val ); }

        vector<PhenotypeNode *>::iterator phenotype_begin() { return individuals->getPhenotypeTree()->node_begin(); }
        vector<PhenotypeNode *>::iterator phenotype_end() { return individuals->getPhenotypeTree()->node_end(); }

        // MarkerCollection Wrapper functions
        const Marker *createMarker( string &id, string &chrom, uint start, uint end, double gPos, string &alleles ) { return markers->createMarker( id, chrom, start, end, gPos, alleles); }
        int getMarkerIndex( const Marker * m ) const;

        vector< Marker *>::iterator marker_begin() { return markers->marker_begin(); }
        vector< Marker *>::iterator marker_end() { return markers->marker_end(); }

        void createGenotypedIndividuals( const vector<int> &indices );
        void createGenotypedMarkers( const vector<int> &indices );

        void dump_data();

        void updateGenotypeTable();

        GenoTable * getGenotypeTable() { return geno_tbl; }

        int getGenotypedIndividualIndex( const std::string & id ) const;
        int getGenotypedMarkerIndex( const std::string & id ) const;

        string getGenotypedMarkerID( int order ) { return genotyped_markers->getIDAtOrderedIndex( order ); }
        string getGenotypedIndividualID( int order) { return genotyped_individs->getIDAtOrderedIndex(order); }

        void getRandomGenotypedIndividualIDSet( std::set<std::string> & id_set, int count = 10) { RandomIDSet( genotyped_individs, id_set, count); }
        void getRandomGenotypedMarkerIDSet( std::set<std::string> & id_set, int count = 100 ) { RandomIDSet( genotyped_markers, id_set, count); }

        void setCaseControlSet( set<string> * case_set, set<string> * ctrl_set );
        CaseControlSet * getCaseControlSet() { return ccs; }

        virtual ~GeneticData();
    protected:

    private:
        IndividualCollection *individuals;
        MarkerCollection *markers;
        CaseControlSet *ccs;

        indexer *genotyped_individs, *genotyped_markers;
        indexer *phenotyped_individs, *phenotyped_traits;

        GenoTable * geno_tbl;

        bool phased;
};


}
}

#endif // GENETICDATA_H
