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
#include "algorithms/genetic_data_func.h"

namespace libgwaspp {
namespace algorithms {

void RandomIDSet( indexer *gd, std::set<std::string> & individ_set, int count ) {
    if( count < gd->included_size() ) {
        const gsl_rng_type *T;
        gsl_rng *r;

        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc( T );

        individ_set.clear();

        double rnd;
        int idx, included = gd->included_size();
        set<int> rand_index;


        while(( int )rand_index.size() < count ) {
            rnd = gsl_rng_uniform( r );
            idx = ( int )( included * rnd ) - 1;

            if( idx < 0 ) { idx = included; }

            rand_index.insert( idx );
        }

        for( set<int>::iterator it = rand_index.begin(), end = rand_index.end(); it != end; ++it ) {
            string tmp = gd->getIDAtOrderedIndex( *it );
            individ_set.insert( tmp );
        }

        gsl_rng_free( r );
        rand_index.clear();

    } else {
        for( int i = 0; i < gd->included_size(); ++i) {
            string tmp = gd->getIDAtOrderedIndex( i );
            individ_set.insert( tmp );
        }
    }
}

}
}
