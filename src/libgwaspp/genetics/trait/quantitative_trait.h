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
#ifndef QUANTITATIVETRAIT_H
#define QUANTITATIVETRAIT_H

#include "genetics/trait/trait.h"

namespace libgwaspp {
namespace genetics {


/**
    A QuantitativeTrait is ...

    @author Patrick Putnam
*/
class QuantitativeTrait : public Trait {
    public:
        QuantitativeTrait( string &_name, string &_note, double _herit, double _mean = 0.0, double _var = 0.0, double _min = 0.0, double _max = 0.0 ) :
            Trait( _name, _note, _herit ), mean( _mean ), var( _var ), min( _min ), max( _max ) {};

        double getMean() { return this->mean; }
        double getVariance() { return this->var; }
        double getMin() { return this->min; }
        double getMax() { return this->max; }

        void setMean( double m ) { this->mean = m; }
        void setVariance( double v ) { this->var = v; }
        void setMax( double m ) { this->max = m; }
        void setMin( double m ) { this->min = m; }

        virtual ~QuantitativeTrait() {}
    protected:
        double mean, var;
        double min, max;
    private:
};

}
}
#endif // QUANTITATIVETRAIT_H
