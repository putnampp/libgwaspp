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
#ifndef SINGLE_MARKER_ANALYZABLE_H
#define SINGLE_MARKER_ANALYZABLE_H

#include "common.h"
#include "genetics/analyzable/case_control_set.h"
#include "genetics/genotype/common_genotype.h"

namespace libgwaspp {
namespace genetics {

class GenotypeDistribution {
public:
    GenotypeDistribution() : current_rIdx( -1 ) {
        reset();
    }

    uint getCurrentIndex() {
        return current_rIdx;
    }
    const frequency_table *getDistribution() {
        return &distribution;
    }
    const header_table *getGenotypes() {
        return &genotypes;
    }

    void setDistribution( const frequency_table &ft ) {
        CopyFrequencyTable(distribution, ft);
    }

    void reset() {
//        distribution.l0 = 0;
//        distribution.l1 = 0;
        ResetFrequencyTable( distribution );

        ResetHeaderTable(genotypes);

        current_rIdx = -1;
    }

    virtual ~GenotypeDistribution() {}

    friend class SingleMarkerAnalyzable;

protected:
    uint current_rIdx;
    frequency_table distribution;
    header_table genotypes;
};

class CaseControlGenotypeDistribution {
public:
    CaseControlGenotypeDistribution() : current_rIdx( -1 ) {
        reset();
    }

    uint getCurrentIndex() {
        return current_rIdx;
    }
    const frequency_table *getCaseDistribution() {
        return &case_dist;
    }
    const frequency_table *getControlDistribution() {
        return &control_dist;
    }
    const header_table *getGenotypes() {
        return &genotypes;
    }

    void setCaseDistribution( const frequency_table &ft) {
        CopyFrequencyTable( case_dist, ft);
    }

    void setControlDistribution( const frequency_table &ft) {
        CopyFrequencyTable( control_dist, ft);
    }

    void reset() {
        case_dist.l0 = 0;
        case_dist.l1 = 0;
        control_dist.l0 = 0;
        control_dist.l1 = 0;

        ResetHeaderTable(genotypes);

        current_rIdx = -1;
    }

    virtual ~CaseControlGenotypeDistribution() {}

    friend class SingleMarkerAnalyzable;

protected:
    uint current_rIdx;
    frequency_table case_dist, control_dist;
    header_table genotypes;
};

class SingleMarkerAnalyzable {
public:
    SingleMarkerAnalyzable() : row_selected( false ), current_dist_rIdx( -1 ) {}
    virtual void selectMarker( uint rIdx ) = 0;
    frequency_table *getGenotypeDistribution( ) {
        return &gt_dist;
    }
    header_table *getGenotypeDistributionHeader( ) {
        return &gt_header;
    }

    virtual void getGenotypeDistribution( uint rIdx, GenotypeDistribution &dist ) = 0;
    virtual void getCaseControlGenotypeDistribution( uint rIdx, CaseControlSet &ccs, CaseControlGenotypeDistribution &ccgd ) = 0;

    virtual void getCaseControlGenotypeDistribution( uint rIdx, CaseControlGenotypeDistribution &ccgd ) = 0;

    virtual ~SingleMarkerAnalyzable() {}
protected:
    void resetSelectedMarker() {
        gt_dist.l0 = 0;
        gt_dist.l1 = 0;
        ResetHeaderTable( gt_header );

        current_dist_rIdx = -1;
        row_selected = false;
    }

    bool row_selected;
    header_table gt_header;
    frequency_table gt_dist;

    uint current_dist_rIdx;

};

}
}


#endif // SINGLE_MARKER_ANALYZABLE_H
