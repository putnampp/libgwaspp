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
#include "genetics/individual/individual.h"

namespace libgwaspp {
namespace genetics {

void Individual::setSex( string &s) {
    if( s[0] == 'F' || s[0] == 'f' || s[0] == '0' ) {
        setSex( FEMALE );
    } else if( s[0] == 'M' || s[0] == 'm' || s[0] == '1') {
        setSex( MALE );
    } else {
        setSex( SEX_UNKNOWN );
    }
}

void Individual::setAge( string &s ) {
    if( s == "NA") {
        this->age = -1.0;
    } else {
        istringstream iss( s );
        iss >> std::dec >> this->age;
    }
}

Individual::~Individual() {
    phenotypes.clear();
}

std::ostream &operator<<( std::ostream &os, const Individual &individ ) {
//    char *sep = "-------------------------------------------------------------\n";
//
//	fprintf(fh, "\n%s", sep);
//	fprintf(fh, "%10s: %-10s\n",  "IND_ID", ind->id);
//	fprintf(fh, "%10s: %-10s\n",  "PED_ID", ind->pedigree_id);
//	fprintf(fh, "%10s: %-10s\n",  "FATHER_ID", ind->father_id);
//	fprintf(fh, "%10s: %-10s\n",  "MOTHER_ID", ind->mother_id);
//	fprintf(fh, "%10s: %c\n",     "SEX", ind->sex);
//	fprintf(fh, "%10s: %-4.1f\n", "AGE", ind->age);
//	fprintf(fh, "%10s: %c\n",     "AFFECTED", ind->affected);
//	fprintf(fh, "%10s: %-10s\n",  "POP_ID", ind->population_id);
//	fprintf(fh, "%10s: %d\n",     "GENOTYPED", ind->genotyped);
//	fprintf(fh, "%10s: %d\n",     "PHENOTYPED", ind->phenotyped);
//	fprintf(fh, "%s", sep);

//  printf("%6d %10s %10s %10s %10s %4c %4.1f %9c %10s %10d %11d\n",

    os.width( 10 );
    os << left << individ.id << individ.pedigree_id << individ.father->id << individ.mother->id;
    os.width( 4 );
    Sex s = individ.getSex();

    if( s == MALE )
        os << left << "M";
    else if( s == FEMALE )
        os << left << "F";
    else
        os << left << "?";

    os.width( 6 );
    os.precision( 1 );
    os << left << fixed << individ.age;

    Affected aff = individ.getAffected();

    os.width( 10 );
    if( aff == AFFECTED )
        os << left << "Affected";
    else if( aff == UNAFFECTED )
        os << left << "Unaffected";
    else
        os << left << "Unknown";

    os << left << individ.population_id;
    os << left << individ.isGenotyped();
    os << left << individ.isPhenotyped();

    return os;
}

}
}
