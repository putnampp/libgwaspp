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
#ifndef HOMOMORPHICMAP_H
#define HOMOMORPHICMAP_H

#include <map>

namespace util {

/**
 *
 * A Homomorphic Map is essentially a bi-directional map between
 * a object L(eft) and some indexing set R(ight)
 *
 * The intended use case is to generate an Alphabet for a categorical
 * data.
 *
 * As an example, assume you have a vector of strings. Within this vector
 * only 3 strings appear {"A Big Data", "Bigger Data", "Catastrophically Huge"}.
 * Obviously, repeatly storing unique instances of these strings is memory inefficient.
 * A simple solution to addressing the memory inefficiency is to store a pointer to each string.
 * In other words, use the memory address of the strings as a homomorphism
 **/

template< class L, typename R = unsigned char, class Compare = less< L >, class Compare2 = less< R > >
class HomomorphicMap
{
    public:
        HomomorphicMap();

        R translateLtoR ( L & left );
        L translateRtoL ( R & right );

        virtual ~HomomorphicMap();
    protected:
    private:
};

}

#endif // HOMOMORPHICMAP_H
