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
#include "genetics/genotype/genotype_record.h"

namespace libgwaspp {
namespace genetics {

GenotypeRecord::GenotypeRecord( map< GenotypeID, EncodedID > &h, vector< EncodedID > &d) {
    //ctor
    initialize( h, d );
}

void GenotypeRecord::initialize( map< GenotypeID, EncodedID > &h, vector< EncodedID > &d ) {
    int header_size = (int) h.size() * sizeof(GenotypeID);
    int data_size = (int) d.size() * sizeof(EncodedID);
    int tot_mem_size = header_size + data_size + 1; // + 1 => to gaurantee NULL terminated

    mem_block = new char[ tot_mem_size ];

    header = reinterpret_cast< GenotypeID * >( mem_block );
    data = reinterpret_cast< EncodedID * >(mem_block + header_size);
    data_end = reinterpret_cast< EncodedID * >( mem_block + tot_mem_size -1 );

    *data_end = 0;   // to gaurantee NULL terminated

    for( map< GenotypeID, EncodedID >::iterator it = h.begin(); it != h.end(); it++ ) {
        header[ it->second ] = it->first;
    }

    int idx = 0;
    for( vector< EncodedID >::iterator it = d.begin(); it != d.end(); it++, idx++) {
        data[idx] = *it;
    }
}

GenotypeRecord::~GenotypeRecord() {
    //dtor
    delete [] mem_block;
}

}
}
