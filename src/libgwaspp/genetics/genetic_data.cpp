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
#include "genetics/genetic_data.h"

namespace libgwaspp {
namespace genetics {

GeneticData::GeneticData( bool _phase ) : ccs(NULL), genotyped_individs(NULL), genotyped_markers(NULL), phenotyped_individs(NULL), phenotyped_traits(NULL), geno_tbl(NULL) {
    //ctor
    individuals = new IndividualCollection();
    markers = new MarkerCollection();

    phased = _phase;
}

int GeneticData::getMarkerIndex( const Marker *m ) const {
    return (*markers)(*m->getIDPtr());
}

void GeneticData::updateGenotypeTable() {
//    gt = new GenotypeTable2( genotyped_markers, genotyped_individs );

#if COMPRESSION_LEVEL == 1
    geno_tbl = new CompressedGenotypeTable( genotyped_markers, genotyped_individs );
#elif COMPRESSION_LEVEL == 2
    geno_tbl = new CompressedGenotypeTable2( genotyped_markers, genotyped_individs );
#elif COMPRESSION_LEVEL == 3
    geno_tbl = new CompressedGenotypeTable3( genotyped_markers, genotyped_individs );
#else
    geno_tbl = new BasicGenotypeTable( genotyped_markers, genotyped_individs );
#endif
}

void GeneticData::addGenotype( int r, int c, const string &s ) {
    geno_tbl->addGenotype( r, c, s );
}

const char *GeneticData::getGenotype( int r, int c ) {
    return geno_tbl->decodeGenotype( (*geno_tbl)(r, c) );
}

void GeneticData::createGenotypedIndividuals( const vector<int> &indices ) {
    if( genotyped_individs != NULL ) { delete genotyped_individs; }

    genotyped_individs = new IndexedSet( individuals, indices );

    cout << "Created IndexedSet of " << genotyped_individs->included_size() << " Individuals" << endl;
}

void GeneticData::createGenotypedMarkers( const vector<int> &indices ) {
    if( genotyped_markers != NULL) { delete genotyped_markers; }

    this->genotyped_markers = new IndexedSet( markers, indices );

    cout << "Created IndexedSet of " << genotyped_markers->included_size() << " Markers" << endl;
}

int GeneticData::getGenotypedIndividualIndex( const std::string &id ) const {
    assert( genotyped_individs != NULL );
    return genotyped_individs->orderOf( id );
}

int GeneticData::getGenotypedMarkerIndex( const std::string &id ) const {
    assert( genotyped_markers != NULL );
    return genotyped_markers->orderOf(id);
}

void GeneticData::setCaseControlSet( set<string> * case_set, set<string> * ctrl_set ) {
    if( ccs == NULL ) {
        ccs = new CaseControlSet( this->genotyped_individs );
    }

    ccs->reset();



    set< int > cases, controls;
    set< int >::iterator dup_it;

    int idx;
    if( case_set != NULL ) {
        for( set<string>::iterator it = case_set->begin(); it != case_set->end(); it++) {
            if( (idx = getGenotypedIndividualIndex( *it ) ) >= 0 )
                cases.insert( idx );
            else {
                cout << "Case Individual not genotyped: " << *it << endl;
            }
        }
    }

    if( ctrl_set != NULL ) {
        for( set<string>::iterator it = ctrl_set->begin(); it != ctrl_set->end(); it++ ) {
            if( (idx = getGenotypedIndividualIndex( *it ) ) >= 0 ) {
                controls.insert( idx );
            } else {
                cout << "Control Individual not genotyped: " << *it << endl;
            }
        }
    }

    for( set<int>::iterator it = cases.begin(); it != cases.end(); it++) {
        if( (dup_it = controls.find( *it )) != controls.end() ) {
            cout << "Index: " << *dup_it << " found in both Case/Controls; Removing from controls" << endl;
            controls.erase( dup_it );
        }
    }

    ccs->setCases( cases );
    ccs->setControls( controls );
}

GeneticData::~GeneticData() {
    //dtor

    delete genotyped_individs;
    delete genotyped_markers;
    delete phenotyped_individs;
    delete phenotyped_traits;

    delete geno_tbl;

    delete individuals;
    delete markers;

    delete GenotypeRecordFactory::instance();


}

}
}
