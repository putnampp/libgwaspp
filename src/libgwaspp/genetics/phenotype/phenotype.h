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
#ifndef PHENOTYPE_H
#define PHENOTYPE_H

#include <iostream>
#include <cmath>

//#include "serializable.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

/**
    Templated class for Phenotypes
    Phenotypes can be defined with a varying set of values
    BasicPhenotypes are associated with strings
    QuantitativePhenotypes defined with some specific numeric (double) value.

    Future types of Phenotypes can be implemented in straightforward manner
 */
template< class T >
class Phenotype {
    public:
        /**
            construct a Phenotype with specified value and a given age
        */
        Phenotype( T &val, float a = NAN ) : value( val ), age( a ) {}

        T getValue() const { return value; }
        float getAge() const { return age; }

        virtual ~Phenotype();
    protected:
        T value;
        float age;
};

/**
    Partial implementation of a phenotype for when the
 */
template < class T >
class Phenotype< T *> {
    public:

        Phenotype( T &val, float a = NAN ) : value( val ), age( a ) {}

        T getValue() const { return *value; }
        float getAge() const { return age; }

        ~Phenotype() { delete value; }
    protected:
        T *value;
        float age;
};

/**
    Phenotypes which are defined with string values (Unknown, Disease, Qualitative)
    often should have a type associated with them to simplify downstream analysis purposes.

    Although, if each set of phenotypes were programmatically separated (ie. given own variables instead of a single pool)
    then type would be implied and any downstream analysis for different types could be simplified
*/
template< >
class Phenotype<string> {
    public:
        Phenotype( string &val, float a = NAN, int t = 0 ) : value( val ), age( a ), type( t ) { }

        string getValue() const { return value; }
        float getAge() const { return age; }
        int getType() const { return type; }

    protected:
        string value;
        float age;
        int type;
};

/**
    Used for: Unknown, qualitative, and disease phenotypes

    Arguably, associated a specific type with each phenotype (ie. defining a 'type' variable)
 */
typedef Phenotype<string> BasicPhenotype;

/**
    Quantitative
*/
typedef Phenotype<double> QuantitativePhenotype;

//typedef streambuilder< BasicPhenotype > PhenotypeBuilder;

}
}
#endif // PHENOTYPE_H
