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
#include "genetics/genotype/genotype_record_factory.h"

namespace libgwaspp {
namespace genetics {

GenotypeRecordFactory::GenotypeRecordFactory() {
    //ctor
    gt_collect = new GenotypeCollection();
}

GenotypeRecordFactory::~GenotypeRecordFactory() {
    //dtor
    delete gt_collect;

    rec_header.clear();
    rec_data.clear();
}

GenotypeRecordFactory *GenotypeRecordFactory::instance() {
    static GenotypeRecordFactory *instance = new GenotypeRecordFactory();
    return instance;
}



bool GenotypeRecordFactory::addRecordData( string &gt ) {
    //GenotypeID id = gt_collect->findOrCreateGenotypeID( gt );

    EncodedID e_id = addGenotypeToHeader( gt_collect->findOrCreateGenotypeID( gt ) );
//    cout << "GT: " << (int)e_id << endl;

    rec_data.push_back( e_id );

    return true;
}

EncodedID GenotypeRecordFactory::addGenotypeToHeader( GenotypeID id ) {
//    map< GenotypeID, EncodedID >::iterator header_it = rec_header.find( id );
//    if( header_it == rec_header.end() ) {
//        // new genotype
//        EncodedID e_id = (EncodedID) rec_header.size();
//        pair< map< GenotypeID, EncodedID >::iterator, bool > success =
//            rec_header.insert( pair< GenotypeID, EncodedID >( id, e_id ) );
//
//        return e_id;
//    }
    if( id == last_gt_id ) {
        return last_enc_id;
    }

    for( map< GenotypeID, EncodedID >::iterator header_it = rec_header.begin(); header_it != rec_header.end(); header_it++ ) {
        if( header_it->first == id ) {
            last_gt_id = id;
            last_enc_id = header_it->second;
            return last_enc_id;
        }
    }

    last_gt_id = id;
    last_enc_id = ( EncodedID ) rec_header.size();
    pair< map< GenotypeID, EncodedID >::iterator, bool > success =
        rec_header.insert( pair< GenotypeID, EncodedID >( id, last_enc_id ) );

    return last_enc_id;
}

GenotypeRecord *GenotypeRecordFactory::finalizeRecord() {
    GenotypeRecord *rec = new GenotypeRecord( rec_header, rec_data );

    rec_header.clear();
    rec_data.clear();

    return rec;
}

}
}
