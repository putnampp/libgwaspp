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
#ifndef PAIRWISE_MARKER_ANALYZABLE_H
#define PAIRWISE_MARKER_ANALYZABLE_H

#include <fstream>

#include "common.h"
#include "genetics/analyzable/case_control_set.h"
#include "genetics/genotype/common_genotype.h"
#include "genetics/genotype/common_genotype_func.h"

#include "genetics/genotype/case_control_selectable.h"

namespace libgwaspp {
namespace genetics {

class ContingencyTable {
    public:
        ContingencyTable() : ma_rIdx( -1 ), mb_rIdx( -1 ) { reset(); }

        uint getMarkerAIndex() { return ma_rIdx; }
        uint getMarkerBIndex() { return mb_rIdx; }

        void setMarkerAIndex( uint ma ) { ma_rIdx = ma; }
        void setMarkerBIndex( uint mb ) { mb_rIdx = mb; }

        const CONTIN_TABLE_T * getContingencyTable() { return &cont; }

        const header_table *getMarkerAHeader() { return &ma_header; }
        const header_table *getMarkerBHeader() { return &mb_header; }


        void setContingency( const CONTIN_TABLE_T & ct ) {
            CopyContingencyTable( cont, ct );
        }

        void setHeaderA( const ushort *headA ) {
            ma_header.xx = headA[0];
            ma_header.aa = headA[1];
            ma_header.ab = headA[2];
            ma_header.bb = headA[3];
        }

        void setHeaderB( const ushort *headB ) {
            mb_header.xx = headB[0];
            mb_header.aa = headB[1];
            mb_header.ab = headB[2];
            mb_header.bb = headB[3];
        }

        void setHeaderA( const header_table & ht ) {
            ma_header.l = ht.l;
        }

        void setHeaderB( const header_table & ht ) {
            mb_header.l = ht.l;
        }

        void reset() {
            ResetHeaderTable(ma_header);
            ResetHeaderTable(mb_header);

            ResetContingencyTable( cont );

            ma_rIdx = -1; mb_rIdx = -1;
        }

        virtual ~ContingencyTable() { }

    protected:
        uint ma_rIdx, mb_rIdx;
        CONTIN_TABLE_T cont;
        header_table ma_header, mb_header;
};

class CaseControlContingencyTable {
    public:
        friend ostream& operator<<( ostream & out, CaseControlContingencyTable & ccct );

        CaseControlContingencyTable() : ma_rIdx( -1 ), mb_rIdx( -1 ) { reset(); }

        uint getMarkerAIndex() { return ma_rIdx; }
        uint getMarkerBIndex() { return mb_rIdx; }

        void setMarkerAIndex( uint rIdx ) { ma_rIdx = rIdx; }
        void setMarkerBIndex( uint rIdx ) { mb_rIdx = rIdx; }

        const header_table *getMarkerAHeader() { return &ma_header; }
        const header_table *getMarkerBHeader() { return &mb_header; }

        const CONTIN_TABLE_T * getCaseContingencyTable() { return &case_contin; }
        const CONTIN_TABLE_T * getControlContingencyTable() { return &control_contin; }


        void setCaseContingency( const CONTIN_TABLE_T & ct ) {
            CopyContingencyTable( case_contin, ct );
        }

        void setControlContingency( const CONTIN_TABLE_T & ct ) {
            CopyContingencyTable( control_contin, ct );
        }

        void updateContingencyTables( const CONTIN_TABLE_T &cs, const CONTIN_TABLE_T &ct ) {
          CopyContingencyTable( case_contin, cs );
          CopyContingencyTable( control_contin, ct );
        }

        void setHeaderA( const ushort *headA ) {
            ma_header.aa = headA[1];
            ma_header.ab = headA[2];
            ma_header.bb = headA[3];
        }

        void setHeaderB( const ushort *headB ) {
            mb_header.aa = headB[1];
            mb_header.ab = headB[2];
            mb_header.bb = headB[3];
        }

        void setHeaderA( const header_table & ht ) {
            ma_header.l = ht.l;
        }

        void setHeaderB( const header_table & ht ) {
            mb_header.l = ht.l;
        }

        void reset() {
            ResetHeaderTable(ma_header);
            ResetHeaderTable(mb_header);

            ResetContingencyTable( case_contin );
            ResetContingencyTable( control_contin );

            ma_rIdx = -1; mb_rIdx = -1;
        }
        virtual ~CaseControlContingencyTable() {}
    protected:
        uint ma_rIdx, mb_rIdx;
        CONTIN_TABLE_T case_contin, control_contin;
        header_table ma_header, mb_header;
};

class PairwiseMarkerAnalyzable : public virtual CaseControlSelectable {
    public:
        PairwiseMarkerAnalyzable() : maIdx( -1 ), mbIdx( -1 ) {}

        virtual void selectMarkerPair( uint rIdx1, uint rIdx2 ) = 0;

        CONTIN_TABLE_T *getContingencyTable() { return &_contingency; }
        header_table *getMarkerAHeader() { return &ma_header; }
        header_table *getMarkerBHeader() { return &mb_header; }

        virtual void getContingencyTable( uint rIdx1, uint rIdx2, ContingencyTable &ct ) = 0;
        virtual void getContingencyTable( uint rIdx1, uint rIdx2, ushort * column_set, ContingencyTable &ct ) = 0;
        virtual void getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlSet &ccs, CaseControlContingencyTable &ccct ) = 0;

        // assumes Case Control Sets have already been selected using selectCaseControl
        virtual void getCaseControlContingencyTable( uint rIdx1, uint rIdx2, CaseControlContingencyTable &ccct ) = 0;
        virtual void getCaseControlContingencyTable( uint rIdx1, uint rIdx2, const marginal_information &m1, const marginal_information &m2, CaseControlContingencyTable & ccct ) = 0;

        virtual ~PairwiseMarkerAnalyzable() {}
    protected:
        void resetSelectedMarkerPair() {
            maIdx = -1; mbIdx = -1;

            ResetContingencyTable( _contingency );
            ResetHeaderTable( ma_header );
            ResetHeaderTable( mb_header );
        }

        uint maIdx, mbIdx;

        CONTIN_TABLE_T _contingency;
        header_table ma_header, mb_header;
};

inline ostream & operator<<(ostream & out, CaseControlContingencyTable & ccct ) {
    out << "Cases" << endl;
    out << ccct.case_contin;
    out << "Controls" << endl;
    out << ccct.control_contin;
    return out;
}

}
}

#endif // PAIRWISE_MARKER_ANALYZABLE_H
