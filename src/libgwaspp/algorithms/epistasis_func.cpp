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
#include "algorithms/epistasis_func.h"

namespace libgwaspp {
namespace algorithms {

void epistasis_all( void *input, void *output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );
    ostream & out = ((output == NULL ) ? cout : *reinterpret_cast< ostream * >(output));

    uint marker_count = inp.gd->getGenotypedMarkersCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    uint total = 0;

    ContingencyTable ct;
    const contingency_table &contingency = *ct.getContingencyTable();

#if DEBUG_LEVEL > 1
    header_table header_a = *ct.getMarkerAHeader();
    header_table header_b = *ct.getMarkerBHeader();
#endif

    INIT_LAPSE_TIME;
    const uint *contin = contingency.contin;

    for( uint i = 0; i < 10; ++i ) {
        RECORD_START;
        for( uint j = i + 1; j < marker_count; ++j ) {
            gt.getContingencyTable( i, j, ct );

            total += *contin++;
            total += *contin++;
            total += *contin++;
            total += *contin++;
            total += *contin++;
            total += *contin++;
            total += *contin++;
            total += *contin++;
            total += *contin++;
            contin = contingency.contin;

#if DEBUG_LEVEL > 1
            out << ( int ) i << " x " << ( int ) j << endl;
            printContingencyTable( contingency, out );
#endif

        }
        RECORD_STOP;

        out << "Single Contingencies " << ( int ) i << " x " << ( int )( marker_count - ( i + 1 ) ) << ": ";
        PRINT_LAPSE(out, "" );
        out << endl;
    }
}

void ContingencyDebug( void * input, void * output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );
    ostream & out = ((output == NULL ) ? cout : *reinterpret_cast< ostream * >(output));

    uint marker_count = inp.gd->getGenotypedMarkersCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    ContingencyTable ct;
    const contingency_table &contingency = *ct.getContingencyTable();

    for( uint i = 0; i < marker_count; ++i ) {
        for( uint j = i + 1; j < marker_count; ++j ) {
            gt.getContingencyTable( i, j, ct );

            out << ( int ) i << " x " << ( int ) j << "\n";
            printContingencyTable( contingency, out );
            out << "\n";
        }
    }
}

void ContingencyPerformance( void * input, void * output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );
    ostream & out = ((output == NULL ) ? cout : *reinterpret_cast< ostream * >(output));

    uint marker_count = inp.gd->getGenotypedMarkersCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    ContingencyTable ct;

    INIT_LAPSE_TIME;

    uint i, j;

    for( i = 0; i < marker_count; ++i ) {
        RECORD_START;
        for( j = i + 1; j < marker_count; ++j ) {
            gt.getContingencyTable( i, j, ct );
        }
        RECORD_STOP;

        out << "Single Contingencies " << ( int ) i << " x " << ( int )( marker_count - ( i + 1 ) ) << ": ";
        PRINT_LAPSE( out, "" );
        out << endl;
    }
}

void epistasis_all2( void *input, void *output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );
    ostream & out = ((output == NULL ) ? cout : *reinterpret_cast< ostream * >(output));

    uint marker_count = inp.gd->getGenotypedMarkersCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    CaseControlSet ccs( gt.getColumnSet() );

    set<int> cases, ctrls;

    for( int i = 0; i < 200; ++i ) {
        cases.insert( i << 1 );
        ctrls.insert( (i << 1) + 1 );
    }

    ccs.setCases( cases );
    ccs.setControls( ctrls );

    double ll=0.0, pval = 0.0;                            // log likelihood

    CaseControlContingencyTable ccct;

    const contingency_table &case_contin = *ccct.getCaseContingencyTable();
    const contingency_table &ctrl_contin = *ccct.getControlContingencyTable();

#if DEBUG_LEVEL > 1
    header_table header_a = *ccct.getMarkerAHeader();
    header_table header_b = *ccct.getMarkerBHeader();
#endif

    INIT_LAPSE_TIME;
    RECORD_START;
    for( uint i = 0; i < marker_count; ++i ) {
        for( uint j = i + 1; j < marker_count; ++j ) {
            gt.getCaseControlContingencyTable( i, j, ccs, ccct );

            ll = pairwise_epi_test( case_contin, ctrl_contin );
            pval = pchisq(ll, 4.0, 0, 0);

#if DEBUG_LEVEL > 1
            out << ( int ) i << " x " << ( int ) j << endl;
            out << "Cases\n";
            printContingencyTable( case_contin, out );

            cout << "\nControls\n";
            printContingencyTable( ctrl_contin, out );
            printf("\nLog likelihood: %f; p-value: %g\n\n", ll, pval);
#endif


        }
    }
    RECORD_STOP;
    out << "Case/Control Contingencies " << ( int )( marker_count ) << ": ";
    PRINT_LAPSE( cout, "" );
    cout << endl;
}

void EpistasisDebug( void * input, void * output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );
    ostream & out = ((output == NULL ) ? cout : *reinterpret_cast< ostream * >(output));

    uint marker_count = inp.gd->getGenotypedMarkersCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    CaseControlSet ccs( gt.getColumnSet() );

    set<int> cases, ctrls;

    for( int i = 0; i < 200; ++i ) {
        cases.insert( i << 1 );
        ctrls.insert( (i << 1) + 1 );
    }

    ccs.setCases( cases );
    ccs.setControls( ctrls );

    double ll=0.0, pval = 0.0;                            // log likelihood

    CaseControlContingencyTable ccct;

    const contingency_table &case_contin = *ccct.getCaseContingencyTable();
    const contingency_table &ctrl_contin = *ccct.getControlContingencyTable();

    for( uint i = 0; i < marker_count; ++i ) {
        for( uint j = i + 1; j < marker_count; ++j ) {
            gt.getCaseControlContingencyTable( i, j, ccs, ccct );

            ll = pairwise_epi_test( case_contin, ctrl_contin );
            pval = pchisq(ll, 4.0, 0, 0);

            out << ( int ) i << " x " << ( int ) j << endl;
            out << "Cases\n";
            printContingencyTable( case_contin, out );

            out << "\nControls\n";
            printContingencyTable( ctrl_contin, out );
            printf("\nLog likelihood: %f; p-value: %g\n\n", ll, pval);
        }
    }
}

void EpistasisPerformance( void * input, void * output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );
    ostream & out = ((output == NULL ) ? cout : *reinterpret_cast< ostream * >(output));

    uint marker_count = inp.gd->getGenotypedMarkersCount();
    GenoTable &gt = *inp.gd->getGenotypeTable();

    CaseControlSet ccs( gt.getColumnSet() );

    set<int> cases, ctrls;

    for( int i = 0; i < 200; ++i ) {
        cases.insert( i << 1 );
        ctrls.insert( (i << 1) + 1 );
    }

    ccs.setCases( cases );
    ccs.setControls( ctrls );

    double ll=0.0, pval = 0.0;                            // log likelihood

    CaseControlContingencyTable ccct;

    const contingency_table &case_contin = *ccct.getCaseContingencyTable();
    const contingency_table &ctrl_contin = *ccct.getControlContingencyTable();

    INIT_LAPSE_TIME;
    RECORD_START;
    for( uint i = 0; i < marker_count; ++i ) {
        for( uint j = i + 1; j < marker_count; ++j ) {
            gt.getCaseControlContingencyTable( i, j, ccs, ccct );

            ll = pairwise_epi_test( case_contin, ctrl_contin );
            pval = pchisq(ll, 4.0, 0, 0);
        }
    }
    RECORD_STOP;
    out << "Case/Control Contingencies " << ( int )( marker_count ) << ": ";
    PRINT_LAPSE(cout, "" );
    cout << endl;
}

void computeBoost( GeneticData * gd, ostream * out ) {
    CaseControlSet &ccs = *gd->getCaseControlSet();

    int nCases = ccs.getCaseCount();
    int nControls = ccs.getControlCount();
    int nIndivids = nCases + nControls;

    GenoTable &gt = *gd->getGenotypeTable();

    // pre-select case/control individuals
    gt.selectCaseControl( ccs );

    // pre-compute marginal distributions for all markers
    marginal_information * pMargins = NULL;
    computeMargins( gt, nIndivids, pMargins );
}

void computeMargins( GenoTable & gt, int nIndivids, marginal_information *& pMargins ) {
    int markerCount = gt.row_size();

    if( pMargins != NULL ) {
        delete [] pMargins;
    }

    pMargins = new marginal_information[ markerCount ];

    CaseControlGenotypeDistribution ccgd;

    for ( int i = 0; i < markerCount; ++i ) {
        // assumes individuals have already been separated into cases and controls
        gt.getCaseControlGenotypeDistribution( i, ccgd, pMargins[i] );
    }
}

double pairwise_epi_test ( const contingency_table &cs, const contingency_table &ct) {
    static int cn[ 9 ];                                 // two-locus genotype count in all samples
    static double pab[ 9 ];                             // conditional genotype probability p(A|B)
    static double pbs[GT_COUNT], pbt[GT_COUNT];                   // conditional genotype probability of the second marker p(B|C)
    static double psa[GT_COUNT], pta[GT_COUNT];                   // conditional case/control probability given the first marker p(C|A)

    int cs1[GT_BUFFER_SIZE], cs2[GT_BUFFER_SIZE];     // genotype count of the two markers in cases
    int ct1[GT_BUFFER_SIZE], ct2[GT_BUFFER_SIZE];     // genotype count of the two markers in controls
    int c1[GT_BUFFER_SIZE], c2[GT_BUFFER_SIZE];       // genotype count of the first and second marker

    int ns=0, nt=0, n=0;                  // total count in cases, controls and all sampels

    double ps = 0, pt = 0, tao=0.0;                   // intermediate variables for likelihood calculation
    double ll=0.0;                            // log likelihood

    int i;                                // genotype of marker 1 and 2

    memset( cs1, 0, GT_COUNT * sizeof(int));
    memset( cs2, 0, GT_COUNT * sizeof(int));
    memset( ct1, 0, GT_COUNT * sizeof(int));
    memset( ct2, 0, GT_COUNT * sizeof(int));
    memset( c1, 0, GT_COUNT * sizeof(int));
    memset( c2, 0, GT_COUNT * sizeof(int));

    // Step 1. calculate marginal counts
    // NOTE: marginal counts for single marker (cs1, cs2, ct1, ct2, c1, c2, ns, nt, n) can be precalculated to increase speed
    int *_cn = cn;
    const uint *_cs = cs.contin, *_ct = ct.contin;
    for( i = 0; i < 9; ) {
        *_cn = *_cs + *_ct;
        switch( i ) {
        case 0:
        case 1:
        case 2:
            cs1[0] += *_cs;
            ct1[0] += *_ct;
            c1[0] += *_cn;
            break;
        case 3:
            ns += cs1[0];
            nt += ct1[0];
            n += c1[0];
        case 4:
        case 5:
            cs1[1] += *_cs;
            ct1[1] += *_ct;
            c1[1] += *_cn;
            break;
        case 6:
            ns += cs1[1];
            nt += ct1[1];
            n += c1[1];
        default:
            cs1[2] += *_cs;
            ct1[2] += *_ct;
            c1[2] += *_cn;
            break;
        }

        switch( i++ ) {
        case 0:
        case 3:
        case 6:
            cs2[0] += *_cs++;
            ct2[0] += *_ct++;
            c2[0] += *_cn++;
            break;
        case 1:
        case 4:
        case 7:
            cs2[1] += *_cs++;
            ct2[1] += *_ct++;
            c2[1] += *_cn++;
            break;
        case 8:
            ns += cs1[2];
            nt += ct1[2];
            n += c1[2];
        default:
            cs2[2] += *_cs++;
            ct2[2] += *_ct++;
            c2[2] += *_cn++;
            break;
        }
    }

    double *_pab = pab;
    _cn = cn;
    for( i = 0; i < 9; ) {
        switch( i ) {
        case 0:
        case 3:
        case 6:
            *_pab = (double) *_cn++ / c2[0];
            break;
        case 1:
        case 4:
        case 7:
            *_pab = (double) *_cn++ / c2[1];
            break;
        default:
            *_pab = (double) *_cn++ / c2[2];
            break;
        }

        switch( i++ ) {
        case 2:
            pbs[0] = (double) cs2[0] / ns;
            pbt[0] = (double) ct2[0] / nt;
            psa[0] = (double) cs1[0] / c1[0];
            pta[0] = 1.0 - psa[0];
            break;
        case 5:
            pbs[1] = (double) cs2[1] / ns;
            pbt[1] = (double) ct2[1] / nt;
            psa[1] = (double) cs1[1] / c1[1];
            pta[1] = 1.0 - psa[1];
            break;
        case 8:
            pbs[2] = (double) cs2[2] / ns;
            pbt[2] = (double) ct2[2] / nt;
            psa[2] = (double) cs1[2] / c1[2];
            pta[2] = 1.0 - psa[2];
            break;
        default:
            break;
        }
    }

    // 3. calculate likelihood (consists of three parts)
    // NOTE: handling of cells with zero counts should be reconsidered
    _cs = cs.contin;
    _ct = ct.contin;
    _pab = pab;

    double *_pbs = pbs, *_pbt = pbt, *_psa = psa, *_pta = pta;
    for( i = 0; i < 9; ++i) {
        if (*_cs > 0) ll += *_cs * log((double) *_cs / n);
        if (*_ct > 0) ll += *_ct * log((double) *_ct / n);

        // part two
        ps = *_pab * *_pbs++ * *_psa;
        pt = *_pab++ * *_pbt++ * *_pta;
        tao += ps + pt;

        switch ( i ) {
        case 2:
        case 5:
            ++_psa;
            ++_pta;
            _pbs = pbs;
            _pbt = pbt;
            break;
        default:
            break;
        }

        if (ps > 0.0)
            ll -= *_cs * log(ps);
        ++_cs;
        if (pt > 0.0)
            ll -= *_ct * log(pt);
        ++_ct;
    }

    // part three
    ll += n*log(tao);

    return 2.0*ll;
}

}
}
