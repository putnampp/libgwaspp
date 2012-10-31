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
 #ifndef PHENOTYPETREE_H
#define PHENOTYPETREE_H

#include <iostream>
#include <set>
#include <vector>
#include <map>

using namespace std;

#include "libgwaspp.h"

namespace libgwaspp {
namespace genetics {

typedef string PhenotypeID;
typedef string PhenotypeValue;

typedef set< const PhenotypeValue *, StringPtrComparer > PhenotypeValues;
typedef set< PhenotypeValue, StringPtrComparer > PhenotypeValuePool;

struct PhenotypeNode {
    PhenotypeID id;
    PhenotypeValues values;

    PhenotypeNode( string &_id ) : id(_id) {}
    ~PhenotypeNode() { values.clear(); }
};

struct PhenotypeNodeComparer {
    bool operator()( const PhenotypeNode & lhs, const PhenotypeNode & rhs ) { return strcmp( lhs.id.c_str(), rhs.id.c_str()) < 0; }
    bool operator()( const PhenotypeNode * lhs, const PhenotypeNode * rhs ) { return (lhs != rhs) && strcmp( lhs->id.c_str(), rhs->id.c_str()) < 0; }
};

class PhenotypeTree {
    public:
        PhenotypeTree();

        PhenotypeNode * addOrGetPhenotype( PhenotypeID &pheno_id );

        const PhenotypeValue * addOrGetPhenotypeValue( PhenotypeID &pheno_id, PhenotypeValue &pheno_val );
        const PhenotypeValue * addOrGetPhenotypeValue( PhenotypeNode * node, PhenotypeValue &pheno_val );

        PhenotypeNode * getPhenotypeAt( int idx ) const { return phenos[ idx ]; }
        int getIndexOf( PhenotypeID &pheno_id );

        int getPhenotypeCount() const { return (int) this->phenos.size(); }

        vector< PhenotypeNode * >::iterator node_begin() { return phenos.begin(); }
        vector< PhenotypeNode * >::iterator node_end() { return phenos.end(); }

        void reset();

        virtual ~PhenotypeTree();
    private:
        typedef map<PhenotypeNode *, int, PhenotypeNodeComparer> PhenoTree;
        PhenoTree lookup_tree;
        vector< PhenotypeNode * > phenos;
        PhenotypeValuePool value_pool;
};

}
}

#endif // PHENOTYPETREE_H
