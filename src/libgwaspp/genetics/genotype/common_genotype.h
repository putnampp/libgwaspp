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
#ifndef COMMON_GENOTYPE_H
#define COMMON_GENOTYPE_H

#include "common.h"

namespace libgwaspp {
namespace genetics {

union genotype_counts {
    uint ui;
    struct {
        byte c0, c1, c2, c3;
    };
};

union header_table {
    ulong l;
    ushort header[ 4 ];
    struct {
        ushort xx, aa, ab, bb;
    };
};

/*
#ifndef ResetHeaderTable
#define ResetHeaderTable(x) (x.l = 0xFFFF000000000000)
#endif
*/

inline void ResetHeaderTable( header_table & ht) {
    ht.l = 0xFFFF000000000000;
}

#define GENOTYPE_COUNT 4 
/// FREQUENCY TABLE
/// Layout:
/// | AA | Aa | aa | xx |
///
union frequency_table {
    uint freq[ GENOTYPE_COUNT ];
    struct {
        ulong l0, l1;
    };
    struct {
        uint aa, ab, bb, xx;
    };
};

inline void CopyFrequencyTable( frequency_table &dest, const frequency_table &src) {
    dest.l0 = src.l0;
    dest.l1 = src.l1;
}

inline void ResetFrequencyTable( frequency_table& ft) {
    ft.l0 = 0;
    ft.l1 = 0;
}

union joint_genotypes {
    char c[10];
    struct {
        byte n0, n1, n2, n3, n4, n5, n6, n7, n8, n9;
    };
    struct {
        byte AA_BB, AA_Bb, AA_bb, Aa_BB, Aa_Bb, Aa_bb, aa_BB, aa_Bb, aa_bb, xx;
    };
    struct {
        ulong ul;
        ushort b;
    };
};

struct marginal_information {
    frequency_table margins, cases, controls;
    double dMarginalEntropy, dMarginalEntropy_Y;
    double dPbc[ 2 * GENOTYPE_COUNT ];   // cases[4]; controls[4]
    double dPca[ 2 * GENOTYPE_COUNT ];   // cases[4]; controls[4]
};

#ifndef ResetJointGenotype
#define ResetJointGenotype( x ) x.ul = 0; x.b = 0;
#endif

#ifndef CopyJointGenotype
#define CopyJointGenotype( dest, src ) dest.ul = src.ul; dest.b = src.b;
#endif

#ifndef SetJointGenotype
#define SetJointGenotype( jg, ma, mb, idx )         \
                    if( ma == 0 || mb == 0 ) {      \
                        ++jg.n9;                    \
                    } else {                        \
                        if( ma == 1 ) {             \
                            if( mb == 1 ) {         \
                                ++jg.n0;            \
                            } else if( mb == 2 ) {  \
                                ++jg.n1;            \
                            } else {                \
                                ++jg.n2;            \
                            }                       \
                        } else if( ma == 2 ) {      \
                            if( mb == 1 ) {         \
                                ++jg.n3;            \
                            } else if( mb == 2 ) {  \
                                ++jg.n4;            \
                            } else {                \
                                ++jg.n5;            \
                            }                       \
                        } else {                    \
                            if( mb == 1 ) {         \
                                ++jg.n6;            \
                            } else if( mb == 2 ) {  \
                                ++jg.n7;            \
                            } else {                \
                                ++jg.n8;            \
                            }                       \
                        }                           \
                    }
#endif

/**********************************************************
 *
 * CONTINGENCY TABLE Definitions
 *
 *********************************************************/

/// Basic Contingency Table
/// Layout:
///   |_BB_|_Bb_|_bb_|
/// AA| 0  | 1  | 2  | 
/// Aa| 3  | 4  | 5  | 
/// aa| 6  | 7  | 8  |
///
/// total xx: 9
union contingency_table {
    uint contin[ 10 ];
    struct {
        uint n0, n1, n2, n3, n4, n5, n6, n7, n8, n9;
    };
    struct {
        uint AA_BB, AA_Bb, AA_bb, Aa_BB, Aa_Bb, Aa_bb, aa_BB, aa_Bb, aa_bb, xx_xx;
    };
};
#define BASIC_CONTIN_COLUMN_COUNT GENOTYPE_COUNT - 1
#define BASIC_CONTIN_ROW_COUNT GENOTYPE_COUNT - 1

/// Extended Contingency Table
/// Layout:
///   |_BB_|_Bb_|_bb_|_xx_|
/// AA| 0  | 1  | 2  | 3  |
/// Aa| 4  | 5  | 6  | 7  |
/// aa| 8  | 9  | 10 | 11 |
/// xx| 12 | 13 | 14 | 15 |
union extended_contingency_table {
    uint contin[ 16 ];
    struct {
        uint AA_BB, AA_Bb, AA_bb, AA_xx, Aa_BB, Aa_Bb, Aa_bb, Aa_xx, aa_BB, aa_Bb, aa_bb, aa_xx,
             xx_BB, xx_Bb, xx_bb, xx_xx;
    };
};
#define EXTENDED_CONTIN_COLUMN_COUNT GENOTYPE_COUNT 
#define EXTENDED_CONTIN_ROW_COUNT GENOTYPE_COUNT

typedef extended_contingency_table CONTIN_TABLE_T;

#ifndef CONTIN_BYTE_SIZE
#define CONTIN_BYTE_SIZE sizeof( CONTIN_TABLE_T )
#endif

#define CONTIN_ROW_COUNT EXTENDED_CONTIN_ROW_COUNT
#define CONTIN_COLUMN_COUNT EXTENDED_CONTIN_COLUMN_COUNT


#ifndef ResetContingencyTable
#define ResetContingencyTable( c ) memset( c.contin, 0, CONTIN_BYTE_SIZE )
#endif

#ifndef CopyContingencyTable
#define CopyContingencyTable( c, v ) memcpy( c.contin, v.contin, CONTIN_BYTE_SIZE ) 
#endif

#ifndef ContingencyAddJointGenotype
#define ContingencyAddJointGenotype( c, v )               \
                        c.AA_BB += v.AA_BB;               \
                        c.AA_Bb += v.AA_Bb;               \
                        c.AA_bb += v.AA_bb;               \
                        c.Aa_BB += v.Aa_BB;               \
                        c.Aa_Bb += v.Aa_Bb;               \
                        c.Aa_bb += v.Aa_bb;               \
                        c.aa_BB += v.aa_BB;               \
                        c.aa_Bb += v.aa_Bb;               \
                        c.aa_bb += v.aa_bb;               \
                        c.xx_xx += v.xx;
#endif

#ifndef UpdateContingency
#define UpdateContingency( c, ma, mb )              \
                    if( ma == 0 || mb == 0 ) {      \
                        ++c.xx_xx;                     \
                    } else {                        \
                        if( ma == 1 ) {             \
                            if( mb == 1 ) {         \
                                ++c.AA_BB;             \
                            } else if( mb == 2 ) {  \
                                ++c.AA_Bb;             \
                            } else {                \
                                ++c.AA_bb;             \
                            }                       \
                        } else if( ma == 2 ) {      \
                            if( mb == 1 ) {         \
                                ++c.Aa_BB;             \
                            } else if( mb == 2 ) {  \
                                ++c.Aa_Bb;             \
                            } else {                \
                                ++c.Aa_bb;             \
                            }                       \
                        } else {                    \
                            if( mb == 1 ) {         \
                                ++c.aa_BB;             \
                            } else if( mb == 2 ) {  \
                                ++c.aa_Bb;             \
                            } else {                \
                                ++c.aa_bb;             \
                            }                       \
                        }                           \
                    }
#endif

#ifndef HeaderStateMachine
#define HeaderStateMachine( geno_code, enc, head_shift, tmp_e, clear_code, head_val ) \
                if( geno_code == 0x0000 ) {                                 \
                    if( isGenotypeHomozygous( enc ) ) {                     \
                        geno_code = 0x1000;                                 \
                        head_shift = 8;                                     \
                        tmp_e = 1;                                          \
                    } else {                                                \
                        geno_code = 0x2000;                                 \
                        head_shift = 4;                                     \
                        tmp_e = 2;                                          \
                    }                                                       \
                    clear_code = 0x0000;                                    \
                } else if( geno_code == 0x1000 ) {                          \
                    if( isGenotypeHomozygous( enc ) ) {                     \
                        geno_code = 0x3000;                                 \
                        head_shift = 0;                                     \
                        tmp_e = 3;                                          \
                    } else {                                                \
                        geno_code = 0x4000;                                 \
                        head_shift = 4;                                     \
                        tmp_e = 2;                                          \
                    }                                                       \
                    clear_code = 0x0F00;                                    \
                } else if( geno_code == 0x2000 ) {                          \
                    assert( isGenotypeHomozygous( enc ) );                  \
                    geno_code = 0x4000;                                     \
                    clear_code = 0x00F0;                                    \
                    head_shift = 8;                                         \
                    tmp_e = 1;                                              \
                } else if( geno_code == 0x3000 ) {                          \
                    assert( !isGenotypeHomozygous( enc ) );                 \
                    geno_code = 0x7000;                                     \
                    clear_code = 0x0F0F;                                    \
                    head_shift = 4;                                         \
                    tmp_e = 2;                                              \
                } else if ( geno_code == 0x4000 ) {                         \
                    assert( isGenotypeHomozygous( enc ) );                  \
                    geno_code = 0x7000;                                     \
                    clear_code = 0x0FF0;                                    \
                    head_shift = 0;                                         \
                    tmp_e = 3;                                              \
                } else {                                                    \
                    assert( false );                                        \
                }                                                           \
                head_val &= clear_code;                                     \
                head_val |= (( geno_code ) | ( enc << head_shift ) );
#endif

#ifndef ParseHeader
#define ParseHeader( _head, _code, _val )                       \
                if( _code == 0x7000 ) {                         \
                    _head.aa = (( _val & 0x0F00 ) >> 8 );       \
                    _head.ab = (( _val & 0x00F0 ) >> 4 );       \
                    _head.bb = ( _val & 0x000F );               \
                } else if( _code == 0x4000 ) {                  \
                    _head.aa = (( _val & 0x0F00 ) >> 8 );       \
                    _head.ab = (( _val & 0x00F0 ) >> 4 );       \
                    _head.bb = 0xFFFF;                          \
                } else if( _code == 0x3000 ) {                  \
                    _head.aa = (( _val & 0x0F00 ) >> 8 );       \
                    _head.ab = 0xFFFF;                          \
                    _head.bb = ( _val & 0x000F );               \
                } else if( _code == 0x2000 ) {                  \
                    _head.aa = 0xFFFF;                          \
                    _head.ab = (( _val & 0x00F0 ) >> 4 );       \
                    _head.bb = 0xFFFF;                          \
                } else if( _code == 0x1000 ) {                  \
                    _head.aa = (( _val & 0x0F00 ) >> 8 );       \
                    _head.ab = 0xFFFF;                          \
                    _head.bb = 0xFFFF;                          \
                } else {                                        \
                    _head.aa = 0xFFFF;                          \
                    _head.ab = 0xFFFF;                          \
                    _head.bb = 0xFFFF;                          \
                }

#endif

}
}

#endif // COMMON_GENOTYPE_H
