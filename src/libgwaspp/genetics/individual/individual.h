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
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <iostream>
#include <iomanip>
#include <map>
#include <cstring>
#include <sstream>

#include "libgwaspp.h"
#include "genetics/phenotype/phenotype_tree.h"

using namespace std;

namespace libgwaspp {
namespace genetics {

enum Phenotyped { UNPHENOTYPED = 0x00, PHENOTYPED = 0x01 };
enum Genotyped { UNGENTOYPED = 0x00, GENOTYPED = 0x02 };
enum Affected { MISSING = 0x00, AFFECTED = 0x10, UNAFFECTED = 0x20 };
enum Sex { SEX_UNKNOWN = 0x00, MALE = 0x40, FEMALE = 0x80 };

typedef vector< const PhenotypeValue * > PhenotypeCollection;

/**
    An Individual is a ...

    @author Patrick Putnam
*/
class Individual {
    public:
        Individual( string _id, PhenotypeTree * pt = NULL ) : id( _id ),  pedigree_id( "" ), population_id( "" ), father( NULL ), mother( NULL ), age( 0. ), ptree(pt), supp_data( 0 ) {}

        string getID() const { return this->id; }
        const string * getIDPtr() const { return &(this->id); }
        string getFatherID() const { return f_id; }
        string getMotherID() const { return m_id; }
        string getPedigreeID() const { return pedigree_id; }
        string getPopulationID() const { return population_id; }
        float getAge() const { return age; }

        void setMother( Individual * m )  { this->mother = m; }
        void setMotherID( const string & mid ) { this->m_id = mid; }

        void setFather( Individual * f )  { this->father = f; }
        void setFatherID( const string & fid ) { this->f_id = fid; }

        void setPedigree( string &ped )  { this->pedigree_id = ped; }
        void setPopulation( string &pop )  { this->population_id = pop; }

        void setAge( float a)  { this->age = a; }
        void setAge( string &s);

//        bool areSiblings( Individual * s ) const { return this->mother == s->mother && this->father == s->father; }

        /**
            Phenotyped information
        */
        inline bool isPhenotyped() const { return ( supp_data & PHENOTYPE ); }
        void setPhenotyped() { unsetPhenotyped(); this->supp_data |= PHENOTYPED; }
        inline void unsetPhenotyped() { this->supp_data &= CLR_PHENO; }

        /**
            Genotyped information
        */
        inline bool isGenotyped() const { return ( supp_data & GENOTYPE ); }
        void setGenotyped() { unsetGenotyped(); this->supp_data |= GENOTYPED; }
        inline void unsetGenotyped() { this->supp_data &= CLR_GENO; }

        /**
            Affectedness information
        */
        Affected getAffected() const { return ( Affected )( supp_data & AFFECTEDNESS ); }
        inline bool isAffected() { return ( supp_data & AFFECTEDNESS ) == AFFECTED; }
        inline bool isUnaffected() { return ( supp_data & AFFECTEDNESS ) == UNAFFECTED; }
        void setAffected( Affected aff ) { clearAffected(); this->supp_data |= ( unsigned char ) aff; }
        inline void clearAffected() { this->supp_data &= CLR_AFFECTED; }

        /**
            Sex information
        */
        Sex getSex() const { return ( Sex )( supp_data & SEX ); }
        void setSex( string &s);
        void setSex( Sex s ) { clearSex(); this->supp_data |= ( unsigned char ) s;}

        inline bool isFemale() { return ( supp_data & SEX ) == FEMALE; }
        inline bool isMale() { return ( supp_data & SEX ) == MALE; }
        inline void clearSex() { this->supp_data &= CLR_SEX; }

        bool operator< (const Individual& b) const { return strcmp(this->id.c_str(), b.id.c_str()) < 0; }

        friend std::ostream &operator<<( std::ostream &os, const Individual &individ );

        const PhenotypeValue * getPhenotypeAt( int idx ) const { return phenotypes[ idx ]; }
        void addPhenotype( const PhenotypeValue *pval ) { this->phenotypes.push_back( pval ); }

        virtual ~Individual();

    protected:
        string id, pedigree_id, population_id;
        string m_id, f_id;
        Individual *father, *mother;

        float age;

        PhenotypeCollection phenotypes;
        PhenotypeTree * ptree;

        unsigned char supp_data;    // SEX | AFFECTEDNESS | GENOTYPED | PHENOTYPED - little-endian bit order [ 7 6 | 5 4 | 3 2 | 1 0 ]
    private:

        enum Mask { PHENOTYPE = 0x01, GENOTYPE = 0x02, AFFECTEDNESS = 0x30, SEX = 0xC0 };
        enum Clear { CLR_PHENO = 0xFE, CLR_GENO = 0xFD, CLR_AFFECTED = 0xCF, CLR_SEX = 0x3F};
};

struct IndividualPtrComparer {
    bool operator() (const Individual * lhs, const Individual * rhs) const { return (*lhs < *rhs); }
};

}
}

#endif // INDIVIDUAL_H
