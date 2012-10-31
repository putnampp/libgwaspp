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

union frequency_table {
    uint freq[ 4 ];
    struct {
        ulong l0, l1;
    };
    struct {
        uint xx, aa, ab, bb;
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
        ulong ul;
        ushort b;
    };
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

union contingency_table {
    uint contin[ 10 ];
    struct {
        ulong l0, l1, l2, l3, l4;
    };
    struct {
        uint n0, n1, n2, n3, n4, n5, n6, n7, n8, n9;
    };
};

#ifndef ResetContingencyTable
#define ResetContingencyTable( c ) c.l0 = 0; c.l1 = 0; c.l2 = 0; c.l3 = 0; c.l4 = 0;
#endif

#ifndef CopyContingencyTable
#define CopyContingencyTable( c, v ) c.l0 = v.l0; c.l1 = v.l1; c.l2 = v.l2; c.l3 = v.l3; c.l4 = v.l4;
#endif

#ifndef ContingencyAddJointGenotype
#define ContingencyAddJointGenotype( c, v )               \
                        c.n0 = c.n0 + v.n0;               \
                        c.n1 = c.n1 + v.n1;               \
                        c.n2 = c.n2 + v.n2;               \
                        c.n3 = c.n3 + v.n3;               \
                        c.n4 = c.n4 + v.n4;               \
                        c.n5 = c.n5 + v.n5;               \
                        c.n6 = c.n6 + v.n6;               \
                        c.n7 = c.n7 + v.n7;               \
                        c.n8 = c.n8 + v.n8;               \
                        c.n9 = c.n9 + v.n9;
#endif

#ifndef UpdateContingency
#define UpdateContingency( c, ma, mb )              \
                    if( ma == 0 || mb == 0 ) {      \
                        ++c.n9;                     \
                    } else {                        \
                        if( ma == 1 ) {             \
                            if( mb == 1 ) {         \
                                ++c.n0;             \
                            } else if( mb == 2 ) {  \
                                ++c.n1;             \
                            } else {                \
                                ++c.n2;             \
                            }                       \
                        } else if( ma == 2 ) {      \
                            if( mb == 1 ) {         \
                                ++c.n3;             \
                            } else if( mb == 2 ) {  \
                                ++c.n4;             \
                            } else {                \
                                ++c.n5;             \
                            }                       \
                        } else {                    \
                            if( mb == 1 ) {         \
                                ++c.n6;             \
                            } else if( mb == 2 ) {  \
                                ++c.n7;             \
                            } else {                \
                                ++c.n8;             \
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
