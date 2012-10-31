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
#ifndef ALLELECOLLECTION_H
#define ALLELECOLLECTION_H

#include <map>
#include <vector>

#include "libgwaspp.h"
#include "genetics/marker/allele_form.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

typedef byte ALLELE_INDEX;

class AlleleCollection {
    public:
        AlleleCollection() {}

        const AlleleForm * findOrCreateAllele( string & alls, char del = 0 );

        ALLELE_INDEX findOrCreateAlleleIndex( string & alls, char del = 0 );
        ALLELE_INDEX findIndexOf( string &alls );

        const AlleleForm * getAlleleAt( int idx ) const { return alleles[idx]; }

        int getAlleleCount() const { return (int)alleles.size(); }

        virtual ~AlleleCollection();
    protected:
    private:
        map< ulong, ALLELE_INDEX > key_lookup;
        vector< AlleleForm * > alleles;
};

}
}

#endif // ALLELECOLLECTION_H
