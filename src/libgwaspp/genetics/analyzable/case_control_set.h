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
 #ifndef CASECONTROLSET_H
#define CASECONTROLSET_H

#include <set>
#include <cstring>
#include <cassert>

#include "libgwaspp.h"

#include "util/index_set/indexer.h"

namespace libgwaspp {
namespace genetics {

using namespace std;
using namespace util;

class CaseControlSet
{
    public:
        CaseControlSet( const indexer * _idx);

        const indexer * getPossibleIndices() { return possible_indices; }

        uint getMaximumIndex() const { return max_index; }
        uint getCaseCount() const { return case_count; }
        uint getControlCount() const { return ctrl_count; }
        uint getTotalCount() const { return total_count; }

        void setCases( const set<int> & case_idx );
        void setControls( const set<int> & ctrl_idx );

        void setAllAsCases();
        void setAllAsControls();

        const ushort * control_begin() { return control_set; }
        const ushort * control_end() { return control_set_end; }

        const ushort * case_begin() { return case_set; }
        const ushort * case_end() { return case_set_end; }

        inline bool isCase( uint idx );
        inline bool isControl( uint idx );

        void reset();

        virtual ~CaseControlSet();
    protected:
        const indexer * possible_indices;
        uint block_count, max_index, byte_count;
        ushort * buffer, * control_set, * case_set;
        ushort * control_set_end, *case_set_end;

        uint case_count, ctrl_count, total_count;

};

}
}

#endif // CASECONTROLSET_H
