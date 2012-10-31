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
#ifndef GENETIC_DATA_FILE_H
#define GENETIC_DATA_FILE_H

#include <sstream>

#include "genetics/genetic_data.h"


namespace libgwaspp {
namespace genetics {

class GeneticDataFile {
public:
    GeneticDataFile() {}

    virtual bool populateGeneticData( string &filename, GeneticData   *gd, char delim = '\t' ) = 0;

    virtual ~GeneticDataFile() {}
protected:
    virtual bool parseHeader( istream *iFile, GeneticData *gd, char delim ) = 0;
    virtual bool parseNextMarkerRecord( istream *iFile, GeneticData *gd, char delim ) = 0;
    virtual bool parseNextGenotypeRecord( istream *iFile, GeneticData *gd, char delim ) = 0;

    virtual bool parseNextRecord( istream * iFile, GeneticData *gd, char delim ) = 0;

    istringstream parser;
    string line, tok;
};

}
}

#endif
