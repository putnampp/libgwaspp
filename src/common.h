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
#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

typedef unsigned char byte;     // 0 - (2^8)  - 1
typedef unsigned short ushort;  // 0 - (2^16) - 1
typedef unsigned int uint;      // 0 - (2^32) - 1
typedef unsigned long ulong;    // 0 - (2^64) - 1

#define bitset_byte_offset(idx) ( idx >> 3 )
#define bitset_bit_offset(idx) (idx & 7)

#define neg_bit_offset(idx) ((~(bitset_bit_offset(idx))) & 0x07)

#define bitset_clear_idx(ptr, idx) ((*(ptr + bitset_byte_offset(idx))) &= ((0xFE << bitset_bit_offset(idx)) | (0x7F >> neg_bit_offset(idx))))

#define bitset_set_idx(ptr, idx) ((*(ptr + bitset_byte_offset(idx))) |= (0x01 << bitset_bit_offset(idx)))

#define bitset_isset(ptr, _idx) (((*(ptr + bitset_byte_offset(_idx))) & (0x01 << bitset_bit_offset(_idx))) != 0)

#define is_bit_set(val, offset) val & (1 << offset)

#define byte_count(bit_count) ((bit_count >> 3) + (((bit_count & 7) > 0) ? 1 : 0))

union byte_block_t {
    byte b;
    char c;
};

union short_memory_block_t {
    ushort us;
    short s;

    // interpret memory block as bytes
    byte bbuf[2];
    char cbuf[2];
    struct {
        byte b0, b1;
    };
    struct {
        char c0, c1;
    };
};

union int_memory_block_t {
    uint ui;
    int i;
    float f;

    // interpret memory block as 2 shorts
    short sbuf[2];
    ushort usbuf[2];

    struct {
        ushort uhi, ulo;
    };

    struct {
        short hi, lo;
    };

    // interpret memory block as 4 chars
    char cbuf[4];
    byte bbuf[4];
    struct {
        byte b0, b1, b2, b3;
    };

    struct {
        char c0, c1, c2, c3;
    };
};

union long_memory_block_t {
    ulong ul;
    long l;
    double d;

    // interpret memory block as 2 integers
    int ibuf[2];
    uint uibuf[2];
    struct {
        uint hi, lo;
    };
    struct {
        float fhi, flo;
    };

    // interpret memory block as 4 shorts
    short sbuf[4];
    ushort usbuf[4];
    struct {
        ushort uhi_hi, uhi_lo, ulo_hi, ulo_lo;
    };
    struct {
        ushort hi_hi, hi_lo, lo_hi, lo_lo;
    };


    // interpret memory block as 8 chars
    char cbuf[8];
    byte bbuf[8];
    struct {
        byte b0, b1, b2, b3, b4, b5, b6, b7;
    };
    struct {
        char c0, c1, c2, c3, c4, c5, c6, c7;
    };
};

union short_memory_block2_t {
    ushort us;
    short s;

    byte bbuf[2];
    char cbuf[2];

    struct {
        byte_block_t hi, lo;
    };
};


union int_memory_block2_t {
    uint ui;
    int i;
    float f;

    short sbuf[2];
    ushort usbuf[2];

    struct {
        short_memory_block2_t hi, lo;
    };
};

// byte b0 = long_memory_block2_t.hi.hi.hi.b; (or long_memory_block2_t.uibuf[0].usbuf[0].bbuf[0];)
// byte b7 = long_memory_block2_t.lo.lo.lo.b;
union long_memory_block2_t {
    ulong ul;
    long l;
    double d;

    int ibuf[2];
    uint uibuf[2];

    struct {
        int_memory_block2_t hi, lo;
    };
};

#endif // COMMON_H_INCLUDED
