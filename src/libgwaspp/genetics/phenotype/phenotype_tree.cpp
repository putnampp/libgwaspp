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
#include "genetics/phenotype/phenotype_tree.h"

namespace libgwaspp {
namespace genetics {

PhenotypeTree::PhenotypeTree() {}

PhenotypeNode *PhenotypeTree::addOrGetPhenotype( PhenotypeID &pheno_id ) {
    PhenotypeNode *node = new PhenotypeNode( pheno_id );
    PhenoTree::iterator pt_iter;

    if(( pt_iter = lookup_tree.find( node ) ) == lookup_tree.end() ) {
        // new phenotype
        phenos.push_back( node );
        lookup_tree.insert( pair< PhenotypeNode *, int> ( node, ( int ) phenos.size() ) );
        return node;
    }

    delete node;

    return pt_iter->first;
}

const PhenotypeValue *PhenotypeTree::addOrGetPhenotypeValue( PhenotypeID &pheno_id, PhenotypeValue &pheno_val ) {
    PhenotypeNode *node = new PhenotypeNode( pheno_id );
    PhenoTree::iterator ptr_iter;

    pair< PhenotypeValuePool::iterator, bool> success = value_pool.insert( pheno_val );

    if(( ptr_iter = lookup_tree.find( node ) ) == lookup_tree.end() ) {
        // new node
        phenos.push_back( node );
        lookup_tree.insert( pair< PhenotypeNode *, int> ( node, ( int ) phenos.size() ) );

        node->values.insert( &*( success.first ) );
    } else {
        delete node;    // clean up unused node

        ptr_iter->first->values.insert( &*( success.first ) );
    }
    return &*( success.first );
}

const PhenotypeValue *PhenotypeTree::addOrGetPhenotypeValue( PhenotypeNode *node, PhenotypeValue &pheno_val ) {
    pair< PhenotypeValuePool::iterator, bool> success = value_pool.insert( pheno_val );
    node->values.insert( &*( success.first ) );
    return &*( success.first );
}

int PhenotypeTree::getIndexOf( PhenotypeID &pheno_id ) {
    PhenotypeNode node( pheno_id );
    PhenoTree::iterator pt_iter;
    if(( pt_iter = lookup_tree.find( &node ) ) == lookup_tree.end() ) {
        return -1;
    }
    return pt_iter->second;
}

void PhenotypeTree::reset() {
    for( vector< PhenotypeNode * >::iterator iter = phenos.begin(); iter != phenos.end(); iter++ ) {
        delete *iter;
    }
    phenos.clear();
    lookup_tree.clear();
    value_pool.clear();
}

PhenotypeTree::~PhenotypeTree() {
    reset();
}

}
}
