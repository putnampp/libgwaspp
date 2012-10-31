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
#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <iostream>

#include "libgwaspp.h"

namespace libgwaspp {
namespace genetics {

#ifndef MAX_ALLELE_COUNT
#define MAX_ALLELE_COUNT 2
#endif

// GENOTYPE_FORM_COUNT defines the maximum number of possible genotypes for a single (marker, individual), including "UNKNOWN"
#ifndef GENOTYPE_FORM_COUNT

#define GENOTYPE_FORM_COUNT 4
#define BITS_PER_GT 2
#define BLOCK_BYTE_SIZE 1
#define GT_PER_BLOCK 4

#define gt_bit_offset(x) ((x << 1) & 7)
#define gt_byte_offset(x) (x >> 2)

#define gt_unmask(val, bit_offset) ((val >> bit_offset) & 3)
#define get_gt_value_at( row_ptr, col ) gt_unmask((*( (row_ptr) + gt_byte_offset(col) )), gt_bit_offset(col))

/*
#define clear_gt_value_at( row_ptr, col) { \
                                            byte bit_offset = gt_bit_offset( col ); \
                                            byte mask = (0xFC << bit_offset) | (0x3F >> (7 - bit_offset)); \
                                            *( (row_ptr) + gt_byte_offset(col) ) &= mask; \
                                         }
*/

#define clear_gt_value_at( row_ptr, col) (*( (row_ptr) + gt_byte_offset(col) )) &= (0xFC << gt_bit_offset( col )) | (0x3F >> (6 - gt_bit_offset( col )))

#define set_gt_value_at( row_ptr, col, val )  { \
                                                byte * ptr = (row_ptr); \
                                                byte bit_offset = gt_bit_offset( col ); \
                                                byte mask = (0xFC << bit_offset) | (0x3F >> (6 - bit_offset)); \
                                                ptr += gt_byte_offset(col); \
                                                (*ptr) &= mask; \
                                                (*ptr) |= (val << bit_offset); \
                                              }

#else
#error Incomplete implementation of user-defined Genotype records
#endif

#ifndef GENOTYPE_ALPHABET
#define GENOTYPE_ALPHABET "ACGT"
#endif

//typedef std::string Genotype;
typedef byte EncodedID;
typedef byte GenotypeID;

enum GENOTYPE_TYPE { GT_UNKNOWN = 0, NA, HETEROZYGOUS, HOMOZYGOUS };

struct Genotype {
    std::string geno;
    GENOTYPE_TYPE type;

    Genotype( std::string & g, GENOTYPE_TYPE t = GT_UNKNOWN ) : geno(g), type(t) {}
};

struct SNPGenotype : public Genotype {
    SNPGenotype( std::string &g, GENOTYPE_TYPE t = GT_UNKNOWN ) : Genotype( g, t ) {
        if( g.length() == 2) {
            if( g == "00" ) {
                type = NA;
            } else if( g[0] == g[1] ) {
                type = HOMOZYGOUS;
            } else {
                type = HETEROZYGOUS;
            }
        }
    }
};

}
}

#endif // GENOTYPE_H
