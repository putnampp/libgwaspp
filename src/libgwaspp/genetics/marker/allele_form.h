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
#ifndef ALLELEFORM_H
#define ALLELEFORM_H

#include <sstream>
#include <vector>

#include "libgwaspp.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

/**
    NOTE:
        - How should duplicate alleles (CC) be handled? Are they allowed? If not permitted, should exceptions be thrown?
        - Should AlleleForms be order independent? For example, is AG the same as GA?
        - Should there be a mechanism for defining one allele as the "dominant" or "recessive"?
        - Should character case be considered? Or, can just force upper case on all allele forms?
*/

class AlleleForm {
    public:
        /**
            Create AlleleForm from a string of alleles separated by some delimiter
            Example:
                AlleleForm( "AC", 0) will separate A and C alleles
                AlleleForm( "A;C", ';') will separate A and C alleles
        **/
        AlleleForm( string &alleles, char delim = 0 );

        /**
            Very simple function for generating a unique hash id from a string of alleles

            @TODO Need to upgrade this hash function to be more robust.
                  Works well for SNPs, but will likely fail for alleles of longer than a single base
        */
        static ulong computeHashID( string &id );

        ulong getID() const { return id; }

        const string *getAlleleAt( int idx ) const { return &alleles[idx]; }
        int getAlleleCount() const { return ( int ) alleles.size(); }

        void addAllele( string &a );

        vector< string >::const_iterator begin() { return alleles.begin(); }
        vector< string >::const_iterator end() { return alleles.end(); }

        /**
            Assumes allele order and character case dependent equality.
        **/
        bool operator==( const AlleleForm &af );

        virtual ~AlleleForm();
    protected:
        void initialize( string &alleles, char delim );
    private:
        ulong id;
        vector< string > alleles;
};

}
}

#endif // ALLELEFORM_H
