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
#ifndef GENOTYPERECORDFACTORY_H
#define GENOTYPERECORDFACTORY_H

#include <iostream>
#include <map>
#include <vector>

#include "libgwaspp.h"
#include "genetics/genotype/genotype.h"
#include "genetics/genotype/genotype_record.h"
#include "genetics/genotype/genotype_collection.h"

using namespace std;

namespace libgwaspp {
namespace genetics {



class GenotypeRecordFactory {
    public:
        typedef map< const string *, int, StringPtrComparer > LookupTable;

        bool addRecordData( string &gt );
        GenotypeRecord * finalizeRecord();

        const Genotype * decodeGenotype( GenotypeID gt ) { return gt_collect->findGenotype(gt); }

        static GenotypeRecordFactory *instance();
        virtual ~GenotypeRecordFactory();
    protected:
        virtual EncodedID addGenotypeToHeader( GenotypeID id );
    private:
        GenotypeRecordFactory();

        GenotypeCollection *gt_collect;

        map< GenotypeID, EncodedID > rec_header;
        vector< EncodedID > rec_data;

        GenotypeID last_gt_id;
        EncodedID last_enc_id;
};

}
}

#endif // GENOTYPERECORDFACTORY_H
