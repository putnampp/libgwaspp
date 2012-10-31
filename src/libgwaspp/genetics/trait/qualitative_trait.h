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
#ifndef QUALITATIVETRAIT_H
#define QUALITATIVETRAIT_H

#include <iostream>
#include <set>
#include "libgwaspp.h"
#include "genetics/trait/trait.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

/**
    @file genetics/traits/qualitative_trait.h

    A Category represents a sub-division of the QualitativeTrait

    @author Patrick Putnam
*/

class Category {
    public :
        Category( string &nm, double f ) : name( nm ), freq( f ) {}

        string getName() const { return name; }
        double getFrequency() const { return freq; }

        bool operator==( const Category &c ) const { return equals( c.name ); }
        bool equals( const string &n ) const { return name.compare( n ) == 0; }

        bool operator<( const Category &c ) const { return name.compare( c.name ) < 0; }


    protected:
        string name;
        double freq;
};

/**
    @file genetics/traits/qualitative_trait.h

    A Qualitative trait is ....

    @author Patrick Putnam
*/

class QualitativeTrait : public Trait {
    public:
        QualitativeTrait( string &name, string &note, double herit ) : Trait( name, note, herit ), missingIdx( 0 ) {}

        void addCategory( Category &c, bool isMissing = false );
        void addCategory( string &name, double freq, bool isMissing = false);

        void setMissingCategory( byte miss ) { this->missingIdx = miss; }

        void setMissingCategory( string n ) ;
        void setMissingCategory( Category *c ) { this->setMissingCategory( c->getName() ); }


        /**
            @param  cat_name    Category name to be located
            @return -1.0 if category name is not found; otherwise, frequency associated with category
        */
        double getFrequency( string &cat_name );

        std::set< Category >::iterator begin() { return cat.begin(); }
        std::set< Category >::iterator end() { return cat.end(); }

        virtual ~QualitativeTrait();
    protected:
        byte missingIdx;
        std::set< Category > cat;
    private:
};

}
}
#endif // QUALITATIVETRAIT_H
