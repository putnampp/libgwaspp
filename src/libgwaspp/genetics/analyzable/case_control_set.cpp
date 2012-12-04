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
#include "genetics/analyzable/case_control_set.h"

namespace libgwaspp {
namespace genetics {

CaseControlSet::CaseControlSet( const indexer *idx ) : possible_indices( idx ) {
    //ctor
    max_index = possible_indices->included_size() - 1;
    block_count = ( possible_indices->included_size() >> 3 ) + 1;
    stream_block_count = ( possible_indices->included_size() >> 4 ) + 1;

    int block_bit_width = (sizeof(ushort) << 3);
    assert( (PROCESSOR_WORD_SIZE % block_bit_width) == 0 );

    int block_per_pword = (PROCESSOR_WORD_SIZE / block_bit_width);
    if( block_per_pword > 1) {
        if( block_count % block_per_pword ) {
            block_count += ( block_per_pword - (block_count % block_per_pword));
        }

        if( stream_block_count % block_per_pword ) {
            stream_block_count += (block_per_pword - (stream_block_count % block_per_pword));
        }
    }

    byte_count = block_count * sizeof(ushort);
    stream_byte_count = stream_block_count * sizeof(ushort);

    cout << "Maximum possible Case/Controls: " << possible_indices->maximum_size() << endl;

    buffer = new ushort[ block_count << 1 ];
    stream_buffer = new ushort[ stream_byte_count << 1 ];

    control_set = buffer;
    case_set = buffer + block_count;

    stream_control_set = stream_buffer;
    stream_case_set = stream_control_set + stream_byte_count;

    control_set_end = case_set;
    case_set_end = case_set + block_count + 1;

    stream_ctrl_end = stream_case_set;
    stream_cs_end = stream_case_set + stream_byte_count + 1;

    case_count = 0; ctrl_count = 0; total_count = 0;
}

void CaseControlSet::setCases( const set<int> & case_idx ) {
    memset( case_set, 0, byte_count );
    memset( stream_case_set, 0, stream_byte_count );

    set<int>::const_iterator it = case_idx.begin(), it_e = case_idx.end();
    uint block_offset, bit_offset;
    uint s_block_offset, s_bit_offset;

    case_count = 0;
    for( ; it != it_e; ++it ) {
        assert( *it <= ( int )max_index );
        block_offset = ( *it >> 3 );
        bit_offset = (( *it & 7 ) << 1 );

        s_block_offset = ( *it >> 4 );
        s_bit_offset = (*it & 0x0F );

        case_set[ block_offset ] = ( case_set[ block_offset ] | ( 0x0003 << bit_offset ) );
        stream_case_set[ s_block_offset ] = ( stream_case_set[ s_block_offset ] | ( 1 << s_bit_offset ) );
        ++case_count;
    }
    assert(case_count == (uint) case_idx.size());
    case_count = ( uint ) case_idx.size();
    total_count = case_count + ctrl_count;

/*
    cout << "Cases: ";
    for( uint i = 0; i < block_count; ++i) {
        cout << hex << (int)case_set[ i ] << " ";
    }
    cout << endl;
    for( uint i = 0; i < stream_block_count; ++i ) {
        cout << hex << (int)stream_case_set[ i ] << " ";
    }
    cout << dec << endl;
*/
}

void CaseControlSet::setControls( const set<int> & ctrl_idx ) {
    memset( control_set, 0, byte_count );
    memset( stream_control_set, 0, stream_byte_count );

    set<int>::const_iterator it = ctrl_idx.begin(), it_e = ctrl_idx.end();
    uint block_offset, bit_offset;
    uint s_block_offset, s_bit_offset;

    ctrl_count = 0;
    for( ; it != it_e; ++it ) {
        assert( *it <= ( int )max_index );
        block_offset = ( *it >> 3 );
        bit_offset = (( *it & 7 ) << 1 );

        s_block_offset = ( *it >> 4 );
        s_bit_offset = ( *it & 0x0F );

        control_set[ block_offset ] = ( control_set[ block_offset ] | ( 0x0003 << bit_offset ) );
        stream_control_set[ s_block_offset ] = ( stream_control_set[ s_block_offset ] | ( 1 << s_bit_offset ));
        ctrl_count++;
    }
    assert(ctrl_count == (uint) ctrl_idx.size());
    ctrl_count = ( uint ) ctrl_idx.size();
    total_count = case_count + ctrl_count;
/*
    cout << "Controls: ";
    for( uint i = 0; i < block_count; ++i ) {
        cout << hex << (int) control_set[ i ] << " ";
    }
    cout << endl;
    for( uint i = 0; i < block_count; ++i ) {
        cout << hex << (int) stream_control_set[ i ] << " ";
    }
    cout << dec << endl;
*/
}

void CaseControlSet::setAllAsCases() {
    memset( case_set, 0xFF, byte_count );
    memset( control_set, 0x00, byte_count );

    memset( stream_case_set, 0xFF, stream_byte_count);
    memset( stream_control_set, 0x00, stream_byte_count);

    case_count = max_index + 1;
    ctrl_count = 0;
    total_count = case_count;

    case_set[ block_count - 1 ] &=  ( 0xFFFF >> (( (7 - ( max_index & 7 )) << 1 ) ) );

/*
    cout << "Cases: ";
    for( uint i = 0; i < block_count; ++i) {
        cout << hex << (int)case_set[ i ] << " ";
    }
    cout << dec << endl;
*/
}

void CaseControlSet::setAllAsControls() {
    memset( case_set, 0x00, byte_count );
    memset( control_set, 0xFF, byte_count );

    memset( stream_case_set, 0x00, stream_byte_count );
    memset( stream_control_set, 0xFF, stream_byte_count );
    
    ctrl_count = max_index + 1;
    case_count = 0;
    total_count = ctrl_count;

    control_set[ block_count - 1 ] &=  ( 0xFFFF >> (( (7 - ( max_index & 7 )) << 1 ) ) );

    cout << "Controls: ";
    for( uint i = 0; i < block_count; ++i ) {
        cout << hex << (int) control_set[ i ] << " ";
    }
    cout << dec << endl;
}

void CaseControlSet::reset() {
    memset( case_set, 0x00, byte_count );
    memset( control_set, 0x00, byte_count );
    memset( stream_case_set, 0x00, stream_byte_count );
    memset( stream_control_set, 0x00, stream_byte_count );
}

bool CaseControlSet::isCase( uint idx ) {
    return (( idx <= max_index ) && ( case_set[( idx >> 3 )] & ( 0x0003 << (( idx & 7 ) << 1 ) ) ) );
}

bool CaseControlSet::isControl( uint idx ) {
    return (( idx <= max_index ) && ( control_set[( idx >> 3 )] & ( 0x0003 << (( idx & 7 ) << 1 ) ) ) );
}

CaseControlSet::~CaseControlSet() {
    //dtor
    delete [] buffer;
    delete [] stream_buffer;
}

}
}
