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
#ifndef GENOTYPETABLE_H
#define GENOTYPETABLE_H

#include <iostream>
#include <set>
#include <vector>

#include "libgwaspp.h"
#include "genetics/genotype/genotype.h"
#include "genetics/genotype/genotype_exceptions.h"
#include "util/exceptions/exceptions.h"
#include "util/table/table.h"

using namespace std;
using namespace util;

namespace libgwaspp {
namespace genetics {

typedef pair<Genotype, byte> EncodedGenotype;

struct gtTableComp {
    bool operator() (const EncodedGenotype &lhs, const EncodedGenotype &rhs) const {
        return lhs.first.compare(rhs.first) < 0;
    }
};

typedef set< EncodedGenotype, gtTableComp > GenotypeSet;

class GenotypeTable {
    public:
        static GenotypeTable *getInstance();

        void loadGenotypeTable( string &gt_file );

        byte addGenotype( EncodedGenotype &eg) throw();

        int getGenotypeCount() const { return (int)genos.size(); }

        byte findEncodedGenotype( Genotype &gt ) throw();
        Genotype getGenotype( byte encoding ) throw();

        virtual ~GenotypeTable();
    protected:


        void reset();
    private:
        GenotypeTable();

        static GenotypeTable *instance;
        GenotypeSet::iterator gtIterator;

        GenotypeSet genos;
        vector< GenotypeSet::iterator > index;
};

class GenotypeTable2 : public Table< byte > {


};

class GenotypeTableBuilder {
    public:
        static GenotypeTableBuilder(
        virtual GenotypeTableBuilder() {}
    private:
        GenotypeTableBuilder();
};

}
}

#endif // GENOTYPETABLE_H
