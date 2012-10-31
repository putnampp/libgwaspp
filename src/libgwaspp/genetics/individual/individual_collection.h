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
#ifndef INDIVIDUALCOLLECTION_H
#define INDIVIDUALCOLLECTION_H

#include <map>
#include <set>
#include <vector>
#include <boost/unordered_map.hpp>

#include "libgwaspp.h"
#include "genetics/individual/individual.h"
#include "genetics/phenotype/phenotype_tree.h"
#include "util/index_set/indexable.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

class IndividualCollection : public util::indexable {
    public:
        typedef boost::unordered_map< string, int > IndLookupSet;

        IndividualCollection() : ptree( new PhenotypeTree() ) {}

        Individual * findOrCreateIndividual( string & id );
        Individual * getIndividualAt( int idx );

        int operator()( const string & id );
        string operator()( int idx ) const;

        int size() const { return (int) individs.size(); }

        int getCount() const { return (int) individs.size(); }
        int getPhenotypeCount() const { return ptree->getPhenotypeCount(); }

        PhenotypeTree * getPhenotypeTree() const { return ptree; }

        vector< Individual * >::iterator begin();
        vector< Individual * >::iterator end();

        virtual ~IndividualCollection();
    protected:
        void reset();

    private:
        IndLookupSet individ_lookup;

        vector< Individual * > individs;
        PhenotypeTree *ptree;
};

}
}

#endif // INDIVIDUALCOLLECTION_H
