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
#ifndef INDIVIDUALGENOTYPEFILE_H
#define INDIVIDUALGENOTYPEFILE_H

#include <fstream>

#include <vector>
#include <boost/algorithm/string.hpp>

#include "libgwaspp.h"
#include "genetics/genetic_data.h"
#include "genetics/genetic_data_file.h"
#include "genetics/individual/individual_collection.h"
#include "gzstream/gzstream.h"
#include "genetics/genotype/genotype_table2.h"

#include "util/time/timing.h"

using namespace std;
using namespace util;
using namespace boost::algorithm;

namespace libgwaspp {
namespace genetics {

class IndividualGenotypeFile : public GeneticDataFile {
    public:
        IndividualGenotypeFile( );

        bool populateGeneticData( string &filename, GeneticData   *gd, char delim );

        virtual ~IndividualGenotypeFile();
    protected:

        virtual void guessExpectedSizes( istream *iFile, int &mcnt, int &gt_width, char delim );

        virtual bool parseHeader( istream *iFile, GeneticData *gd, char delim );
        virtual bool parseNextMarkerRecord( istream *iFile, GeneticData *gd, char delim );
        virtual bool parseNextGenotypeRecord( istream *iFile, GeneticData *gd, char delim );

        virtual bool parseNextRecord( istream * iFile, GeneticData *gd, char delim ) {return true; }

        int row_idx, col_idx;
        vector<int> individ_indexes, marker_indexes;
};

}
}

#endif // INDIVIDUALGENOTYPEFILE_H
