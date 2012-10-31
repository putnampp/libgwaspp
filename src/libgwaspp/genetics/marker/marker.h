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
#ifndef MARKER_H
#define MARKER_H

#include <iostream>
#include <iterator>
#include <set>

#include "libgwaspp.h"
#include "genetics/chromosome/chromosomal_interval.h"
#include "genetics/marker/allele_form.h"
#include "genetics/marker/allele_collection.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

enum MarkerType { UNKNOWN = 0, SNP, SNV, CNV, INDEL, STR, GENE };

struct MarkerNameComparer;
struct MarkerPositionComparer;

/**
    A Marker is a named chromosomal interval

    @author Patrick Putnam
**/
class Marker : public ChromosomalInterval {
    public:
        Marker(string & n, ChromosomeID cIdx, uint s, uint e, double genPos, const AlleleForm * af) : ChromosomalInterval(cIdx, s, e, genPos), id(n), allele(af)  {}

        string getID() const { return id; }
        const string * getIDPtr() const { return &id; }

        const AlleleForm * getAlleles() const { return allele; }

        friend struct MarkerIDComparer;
        friend struct MarkerPositionComparer;

        virtual ~Marker() {}
    protected:
    private:
        string id;
//        ALLELE_INDEX allele_idx;
        const AlleleForm * allele;
};

struct MarkerIDComparer {
    bool operator()( const Marker & lhs, const Marker & rhs ) { return strcmp( lhs.id.c_str(), rhs.id.c_str()) < 0; }
};

struct MarkerPositionComparer {
    bool operator()( const Marker & lhs, const Marker & rhs ) { return (lhs.chromIdx < rhs.chromIdx) || (lhs.chromIdx == rhs.chromIdx && lhs.getStart() < rhs.getStart()); }

};

}
}
#endif // MARKER_H
