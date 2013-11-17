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
    const CONTIN_TABLE_T &contingency = *ct.getContingencyTable();

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
    const CONTIN_TABLE_T &contingency = *ct.getContingencyTable();

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

    uint i, j, k = 0;

    for( i = 0; i < marker_count; ++i ) {
        //RECORD_START;
        for( j = i + 1; j < marker_count; ++j, ++k ) {
            RECORD_START;
            gt.getContingencyTable( i, j, ct );
            RECORD_STOP;

            out << (int)k;
            PRINT_LAPSE( out, "\t");
            out << endl;
        }
        //RECORD_STOP;

        //out << "Single Contingencies " << ( int ) i << " x " << ( int )( marker_count - ( i + 1 ) ) << ": ";
        //PRINT_LAPSE( out, "" );
        //out << endl;
    }
}

void ContingencyCCPerformance( void * input, void * output ) {
    IndexedInput &inp = *reinterpret_cast<IndexedInput *>( input );
    ostream & out = ((output == NULL ) ? cout : *reinterpret_cast< ostream * >(output));

    GeneticData * gd = inp.gd;
    CaseControlSet &ccs = *gd->getCaseControlSet();

    int nCases = ccs.getCaseCount();
    int nControls = ccs.getControlCount();
    int nIndivids = nCases + nControls;
    int nMarkerCount = 0;

    GenoTable &gt = *gd->getGenotypeTable();

    // pre-select case/control individuals
    gt.selectCaseControl( ccs );

    // pre-compute marginal distributions for all markers
    marginal_information * pMargins = NULL;
    marginal_information * pMar1, * pMar2;
    computeMargins( gt, nIndivids, pMargins, nMarkerCount );

    vector< uint > filteredIndices;

    for( uint i = 0; i < (uint)nMarkerCount; ++i ) {
        filteredIndices.push_back(i);
    }

    vector< uint >::const_iterator itIdx1, itIdx2;
    uint idx, idx2;
    CaseControlContingencyTable ccct;

    const CONTIN_TABLE_T & contin_ca = *ccct.getCaseContingencyTable();
    const CONTIN_TABLE_T & contin_co = *ccct.getControlContingencyTable();

/*
    const uint * pContinCa, * pContinCo, *denom;
    double dPab, *pPbc_ca, *pPbc_co, *pPca_ca, *pPca_co;
    double tao, interMeasure;
    double tmp1, tmp2, tmp3;

    typedef pair< uint, uint > SNPPair;
    typedef pair< SNPPair, double> SNPInteractionPair;
    vector< SNPInteractionPair > passingThreshold;

    double maxInteraction = -99999999, minInteraction = 999999999, thresholdRecord = 30.0;
*/
    unsigned int k = 0;
    INIT_LAPSE_TIME;
    // pre-screening
    for( itIdx1 = filteredIndices.begin(); itIdx1 != filteredIndices.end(); itIdx1++ ) {
        idx = *itIdx1;
        pMar1 = &pMargins[idx];
        for( itIdx2 = itIdx1 + 1; itIdx2 != filteredIndices.end(); itIdx2++, ++k ) {
            idx2 = *itIdx2;
            pMar2 = &pMargins[ idx2 ];
            RECORD_START;
            gt.getCaseControlContingencyTable( idx, idx2, *pMar1, *pMar2, ccct );

            RECORD_STOP;
            out << (int)k;
            PRINT_LAPSE( out, "");
            out << endl;
        }
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

    const CONTIN_TABLE_T &case_contin = *ccct.getCaseContingencyTable();
    const CONTIN_TABLE_T &ctrl_contin = *ccct.getControlContingencyTable();

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

    const CONTIN_TABLE_T &case_contin = *ccct.getCaseContingencyTable();
    const CONTIN_TABLE_T &ctrl_contin = *ccct.getControlContingencyTable();

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

    const CONTIN_TABLE_T &case_contin = *ccct.getCaseContingencyTable();
    const CONTIN_TABLE_T &ctrl_contin = *ccct.getControlContingencyTable();

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
    int nMarkerCount = 0;

    GenoTable &gt = *gd->getGenotypeTable();

    // pre-select case/control individuals
    gt.selectCaseControl( ccs );

    // pre-compute marginal distributions for all markers
    marginal_information * pMargins = NULL;
    marginal_information * pMar1, * pMar2;
    computeMargins( gt, nIndivids, pMargins, nMarkerCount );

    vector< uint > filteredIndices;

    for( uint i = 0; i < (uint)nMarkerCount; ++i ) {
        //if( pMargins[i].margins.xx == 0 )
        filteredIndices.push_back(i);
    }

    vector< uint >::const_iterator itIdx1, itIdx2;
    uint idx, idx2;
    CaseControlContingencyTable ccct;

    const CONTIN_TABLE_T & contin_ca = *ccct.getCaseContingencyTable();
    const CONTIN_TABLE_T & contin_co = *ccct.getControlContingencyTable();

    const uint * pContinCa, * pContinCo, *denom;
    double dPab, *pPbc_ca, *pPbc_co, *pPca_ca, *pPca_co;
    double tao, interMeasure;
    double tmp1, tmp2, tmp3;

    typedef pair< uint, uint > SNPPair;
    typedef pair< SNPPair, double> SNPInteractionPair;
    vector< SNPInteractionPair > passingThreshold;

    double maxInteraction = -99999999, minInteraction = 999999999, thresholdRecord = 30.0;

    *out << "Pre-screening " << filteredIndices.size() << " SNP interactions" << endl;

    INIT_LAPSE_TIME;
    RECORD_START;
    // pre-screening
    for( itIdx1 = filteredIndices.begin(); itIdx1 != filteredIndices.end(); itIdx1++ ) {
        idx = *itIdx1;
    //for( idx = 0; idx < (uint) nMarkerCount; ++idx ) {
        pMar1 = &pMargins[idx];
        for( itIdx2 = itIdx1 + 1; itIdx2 != filteredIndices.end(); itIdx2++ ) {
        //for( int idx2 = idx + 1; idx2 < nMarkerCount; ++idx2 ) {
            idx2 = *itIdx2;
            pMar2 = &pMargins[ idx2 ];
            //RECORD_START;
            gt.getCaseControlContingencyTable( idx, idx2, *pMar1, *pMar2, ccct );

            //RECORD_STOP;
            //PRINT_LAPSE( *out, "");
            //*out << endl;
            
            // less loop overhead
            // no index calculations
            // computes P(A | B)
            //
            //  CONTINGENCY TABLE ORGANIZATION
            //     |__________B___________| MAR1
            //  ___|__AA__|__ AB__|___BB__|
            // | AA|  n0  |   n1  |   n2  |
            //A| AB|  n3  |   n4  |   n5  |
            // | BB|  n6  |   n7  |   n8  |
            // MAR2|      |       |       |
            //
            tao = 0.0; interMeasure = 0.0;
            pContinCa = &contin_ca.contin[0];
            pContinCo = &contin_co.contin[0];
            pPbc_ca = &pMar2->dPbc[0];
            pPbc_co = &pMar2->dPbc[GENOTYPE_COUNT];
            pPca_ca = &pMar1->dPca[0];
            pPca_co = &pMar1->dPca[GENOTYPE_COUNT];
            denom = &pMar2->margins.freq[0];

            //*out << idx << "x" << *itIdx2<<endl;

            for( int i = 0, j = 3; i < 9; ++i ) {
                if( ! j-- ) {
                    j = 2;
                    pPbc_ca = &pMar2->dPbc[0];
                    pPbc_co = &pMar2->dPbc[GENOTYPE_COUNT];
                    pPca_ca++;
                    pPca_co++;
                    denom = &pMar2->margins.freq[0];
                    ++pContinCa;    // skip XX column
                    ++pContinCo;    // skip XX column
                }
                dPab = (double)(*pContinCa + *pContinCo) / (double)*denom++;
                tmp2 = (dPab) * (*pPbc_ca++) * (*pPca_ca);
                tmp3 = (dPab) * (*pPbc_co++) * (*pPca_co);
                tao += tmp2 + tmp3;
                if( *pContinCa > 0 ) {
                    tmp1 = (double) *pContinCa / nIndivids;
                    interMeasure += tmp1 * log(tmp1);

                    if( tmp2 > 0 ) {
                        interMeasure += -tmp1 * log(tmp2);
                    }
                }
                pContinCa++;

                if( *pContinCo > 0 ) {
                    tmp1 = (double) *pContinCo / nIndivids;
                    interMeasure += tmp1 * log(tmp1);
                    if( tmp3 > 0 ) {
                        interMeasure += -tmp1 * log(tmp3);
                    }
                }
                pContinCo++;
            }

            interMeasure = (interMeasure + log(tao)) * nIndivids * 2.0;
            //*out << idx << "x" << *itIdx2 << "\t" << tao << "\t" << interMeasure << endl;


            if( interMeasure > maxInteraction ) {
                maxInteraction = interMeasure;
            }

            if( interMeasure < minInteraction ) {
                minInteraction = interMeasure;
            }

            if( interMeasure > thresholdRecord ) {
                passingThreshold.push_back( SNPInteractionPair( SNPPair( idx, idx2 ), interMeasure) );
            }
        }
    }
    RECORD_STOP;
    PRINT_LAPSE( *out, "");
    *out << endl;

    *out << "Located " << passingThreshold.size() << " potential interactions" << endl;

    *out << "Performing deeper analysis of SNPs" << endl;
    vector< double > zval;
    computeGTest(gt, pMargins, nIndivids, passingThreshold, zval );

    vector< SNPInteractionPair >::const_iterator itPair;
    vector< double >::const_iterator itZ;
    idx = 0;
    for( itPair = passingThreshold.begin(), itZ = zval.begin(); itPair != passingThreshold.end(); itPair++, itZ++ ) {
        if( itPair->second > thresholdRecord ) {
            *out << boost::format( "%7d\t%7d\t%7d\t%f\t%f\t%f\t%f" ) % (idx++) % itPair->first.first % itPair->first.second % 0.0 % 0.0 % itPair->second % *itZ;
            *out << endl;
        }
    } 
}

void computeGTest( GenoTable & gt, marginal_information * pMargins, uint nIndivids, vector< SNPInteractionPair > & passingThreshold, vector< double > & zval ) {

    const int JOINT_GENOTYPE_SIZE = GT_COUNT * GT_COUNT;
    const int MU_SIZE = JOINT_GENOTYPE_SIZE * GT_BUFFER_COUNT; // 3 (genotypes_a) * 3 (genotypes_b) * 2 (case_control)
    double mu[ MU_SIZE ] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    double mu0[ MU_SIZE ] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double mutmp[ MU_SIZE ], mu0tmp[ MU_SIZE ];
    double mu_ij[ JOINT_GENOTYPE_SIZE  ], mu_ik[ GT_BUFFER_SIZE ], mu_jk[ GT_BUFFER_SIZE ];
    double muError, tmp1, tmp2, tmp3, tmp4;

    double * pMu, * pMu0;
    double * pMu_ca, *pMu_co;
    double * pMu0_ca, *pMu0_co;
    double * pMu_ik_ca, * pMu_ik_co;
    double * pMu_jk_ca, * pMu_jk_co;

    CaseControlContingencyTable ccct;

    const CONTIN_TABLE_T & contin_ca = *ccct.getCaseContingencyTable();
    const CONTIN_TABLE_T & contin_co = *ccct.getControlContingencyTable();

    const uint * pContinCa, * pContinCo;
    int idx, idx2;
    double tao, interMeasure;
    marginal_information * pMar1, * pMar2;
    uint AlleleJointDistr[ 8 ];

    // Exact Gtest
    for( vector< SNPInteractionPair >::iterator itPair = passingThreshold.begin();
         itPair != passingThreshold.end(); itPair++ ) {
        idx = itPair->first.first;
        idx2 = itPair->first.second;

        pMar1 = &pMargins[idx];
        pMar2 = &pMargins[idx2];
        gt.getCaseControlContingencyTable( idx, idx2, *pMar1, *pMar2, ccct );

        memcpy( mutmp, mu, MU_SIZE * sizeof(double) );
        memcpy( mu0tmp, mu0, MU_SIZE * sizeof(double) );

        muError = 0.0;
        pMu = mutmp;
        pMu0 = mu0tmp;
        for( int i = 0; i < MU_SIZE; ++i ) {
            muError += abs( *pMu - *pMu0 );
        }

        while( muError > 0.001 ) {
            memcpy( mu0tmp, mutmp, MU_SIZE * sizeof(double));

            pMu = mu_ij;
            pMu_ca = mutmp;
            pMu_co = &mutmp[JOINT_GENOTYPE_SIZE];
            pContinCa = contin_ca.contin;
            pContinCo = contin_co.contin;

            memset( mu_ik, 0, GT_BUFFER_SIZE * sizeof(double) );
            memset( mu_jk, 0, GT_BUFFER_SIZE * sizeof(double) );
            
            pMu_ik_ca = mu_ik;
            pMu_ik_co = &mu_ik[GT_COUNT];
            pMu_jk_ca = mu_jk;
            pMu_jk_co = &mu_jk[GT_COUNT];
            for( int i = 0, j = 3; i < JOINT_GENOTYPE_SIZE; ++i ) {
                if( ! j-- ) {
                    j = 2;
                    pMu_ik_ca++;
                    pMu_ik_co++;
                    pMu_jk_ca = mu_jk;
                    pMu_jk_co = &mu_jk[GT_COUNT];

                    ++pContinCa;
                    ++pContinCo;
                }
                *pMu = (*pMu_ca + *pMu_co);
                if( *pMu > 0 ) {
                        *pMu_ca = *pMu_ca * ( *pContinCa + *pContinCo ) / *pMu;
                        *pMu_co = *pMu_co * ( *pContinCa + *pContinCo ) / *pMu;
                } else {
                    *pMu_ca = 0;
                    *pMu_co = 0;
                }

                *pMu_ik_ca += *pMu_ca;
                *pMu_ik_co += *pMu_co;

                *pMu_jk_ca++ += *pMu_ca++;
                *pMu_jk_co++ += *pMu_co++;

                ++pContinCa;
                ++pContinCo;
            }

            
            pMu_ca = mutmp;
            pMu_co = &mutmp[JOINT_GENOTYPE_SIZE];
            pMu0_ca = mu0tmp;
            pMu0_co = &mu0tmp[JOINT_GENOTYPE_SIZE];
            muError = 0.0;
            //pMu_ik_ca = mu_ik;
            //pMu_ik_co = &mu_ik[GT_COUNT];
            //pMu_jk_ca = mu_jk;
            //pMu_jk_co = &mu_jk[ GT_COUNT ];
            for( int i = 0, j = 0, k = 0, l = 0; i < JOINT_GENOTYPE_SIZE; ++i, ++l ) {
                if ( ! j-- ) {
                    j = 2;
                    if( mu_ik[ k ] > 0 )
                        tmp1 = pMar1->cases.freq[k] / mu_ik[ k ];
                    else
                        tmp1 = 0.0;

                    if( mu_ik[ GT_COUNT + k] > 0 )
                        tmp2 = pMar1->controls.freq[k] / mu_ik[ GT_COUNT + k ];
                    else
                        tmp2 = 0.0;
                    ++k;
                    l = 0;
                }

                if(  mu_jk[l] > 0 )
                    tmp3 = pMar2->cases.freq[l] / mu_jk[l];
                else
                    tmp3 = 0.0;

                if(  mu_jk[GT_COUNT + l] > 0 )
                    tmp4 = pMar2->controls.freq[l] / mu_jk[GT_COUNT + l];
                else
                    tmp4 = 0.0;

                *pMu_ca = *pMu_ca * (tmp1) * tmp3;
                *pMu_co = *pMu_co * (tmp2) * tmp4;

                muError += abs( *pMu_ca++ - *pMu0_ca++ );
                muError += abs( *pMu_co++ - *pMu0_co++ );
            }
        } // end while

        tao = 0.0;
        interMeasure = 0.0;
        pContinCa = contin_ca.contin;
        pContinCo = contin_co.contin;
        pMu_ca = mutmp;
        pMu_co = &mutmp[JOINT_GENOTYPE_SIZE];
        for( int i = 0; i < JOINT_GENOTYPE_SIZE; ++i ) {
            if( i == 3 || i == 6 ) {
                ++pContinCa;
                ++pContinCo;
            }
            if( *pContinCa > 0 ) {
                tmp1 = (double) *pContinCa / nIndivids;
                interMeasure += tmp1 * log(tmp1);
            } else 
                tmp1 = 0.0;

            if( *pMu_ca > 0 ) {
                tmp2 = *pMu_ca / nIndivids;
                interMeasure += -tmp1 * log(tmp2);
                tao += tmp2;
            }

            if( *pContinCo > 0 ) {
                tmp1 = (double) *pContinCo / nIndivids;
                interMeasure += tmp1 * log(tmp1);
            } else
                tmp1 = 0.0;

            if( *pMu_co > 0 ) {
                tmp2 = *pMu_co / nIndivids;
                interMeasure += -tmp1 * log(tmp2);
                tao += tmp2;
            }

            pContinCa++; pContinCo++; pMu_ca++; pMu_co++;
        }

        interMeasure = (interMeasure + log(tao)) * nIndivids * 2.0;
        itPair->second = interMeasure;

        AlleleJointDistr[0] = (contin_ca.contin[0]<<2) + (contin_ca.contin[1]<<1) + (contin_ca.contin[4]<<1) + contin_ca.contin[5];
        AlleleJointDistr[1] = (contin_ca.contin[2]<<2) + (contin_ca.contin[1]<<1) + (contin_ca.contin[6]<<1) + contin_ca.contin[5];
        AlleleJointDistr[2] = (contin_ca.contin[8]<<2) + (contin_ca.contin[9]<<1) + (contin_ca.contin[4]<<1) + contin_ca.contin[5];
        AlleleJointDistr[3] = (contin_ca.contin[10]<<2) + (contin_ca.contin[9]<<1) + (contin_ca.contin[6]<<1) + contin_ca.contin[5];

        AlleleJointDistr[4] = (contin_co.contin[0] << 2) + (contin_co.contin[1] << 1) + (contin_co.contin[4]<< 1) + contin_co.contin[5];
        AlleleJointDistr[5] = (contin_co.contin[2]<<2) + (contin_co.contin[1]<<1) + (contin_co.contin[6]<<1) + contin_co.contin[5];
        AlleleJointDistr[6] = (contin_co.contin[8]<<2) + (contin_co.contin[9]<<1) + (contin_co.contin[4]<<1) + contin_co.contin[5];
        AlleleJointDistr[7] = (contin_co.contin[10]<<2) + (contin_co.contin[9]<<1) + (contin_co.contin[6]<<1) + contin_co.contin[5]; 
    
        double or_aff = log( (double)(AlleleJointDistr[0]*AlleleJointDistr[3])/ (double)(AlleleJointDistr[1]*AlleleJointDistr[2]) );
        double v_aff = 1/(double)AlleleJointDistr[0] + 1/(double)AlleleJointDistr[1] + 1/(double)AlleleJointDistr[2] + 1/(double)AlleleJointDistr[3];
    
        double or_unf = log( (double)(AlleleJointDistr[4]*AlleleJointDistr[7])/ (double)(AlleleJointDistr[5]*AlleleJointDistr[6]) );
        double v_unf = 1/(double)AlleleJointDistr[4] + 1/(double)AlleleJointDistr[5] + 1/(double)AlleleJointDistr[6] + 1/(double)AlleleJointDistr[7];

        zval.push_back( (or_aff - or_unf) / sqrt( v_aff + v_unf) );
    }
}

void computeMargins( GenoTable & gt, int nIndivids, marginal_information *& pMargins, int & nMarkerCount ) {
    nMarkerCount = gt.row_size();

    if( pMargins != NULL ) {
        delete [] pMargins;
    }

    pMargins = new marginal_information[ nMarkerCount ];

    CaseControlGenotypeDistribution ccgd;

    for ( int i = 0; i < nMarkerCount; ++i ) {
        // assumes individuals have already been separated into cases and controls
        gt.getCaseControlGenotypeDistribution( i, ccgd, pMargins[i] );
    }
}

double pairwise_epi_test ( const CONTIN_TABLE_T &cs, const CONTIN_TABLE_T &ct) {
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
