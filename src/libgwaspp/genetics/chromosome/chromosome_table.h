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
#ifndef CHROMOSOMETABLE_H
#define CHROMOSOMETABLE_H

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include "libgwaspp.h"
#include "genetics/chromosome/chromosome.h"
#include "genetics/chromosome/chromosome_exceptions.h"
#include "util/exceptions/exceptions.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

typedef pair<Chromosome, byte> EncodedChromosome;

struct chromComp {
    bool operator()( const EncodedChromosome &lhs, const EncodedChromosome &rhs) {
        return lhs.first < rhs.first;
    }
};

typedef set< EncodedChromosome, chromComp> ChromosomeSet;

/**
    A set of chromosomes

    A Singleton instance is provided for situations where a single
    set of chromosomes is in use.
*/
class ChromosomeTable {
    public:

        static ChromosomeTable *getInstance();

        /**
            Loads a Chromosome Table from specified file.

            Expected file format:
                - TAB-, or comma-delimited file
                - Single line per chromosome
                - zeroth column is index of the chromosome
                - first column is chromosome primary name
                - second column is an integer representing the length (# of bases) of the chromosome
                - (OPTIONAL) third column of pipe '|' seperated aliases for chromosome
            idx{TAB|,}chrom_name{TAB|,}chrom_len({TAB|,}alias|alias)
        */
        void loadChromosomeTable( string &chromTable );

        byte findChromosomeEncoding( string &chrom ) throw();
        byte findChromosomeEncoding( const Chromosome &c ) throw();

        Chromosome getChromosome( string &name);
        Chromosome getChromosome( byte code ) throw();

        virtual ~ChromosomeTable();
    protected:
        byte addChromosome( Chromosome &c ) throw();
        void addChromosomeAlias( Chromosome &c, string &alias) throw();

        void reset();
    private:
        ChromosomeTable();
        static ChromosomeTable *instance;
        //map<byte, Chromosome* > chroms;

        ChromosomeSet chroms;
        ChromosomeSet::iterator chromIter;

        vector< ChromosomeSet::iterator > index;

        map<string, ChromosomeSet::iterator > chrom_aliases;
        map<string, ChromosomeSet::iterator >::iterator aliasIter;
};

}
}
#endif // CHROMOSOMETABLE_H
