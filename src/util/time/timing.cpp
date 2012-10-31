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
#include "util/time/timing.h"

namespace util {

int diff_TIME( TIME &res, TIME &x, TIME &y ) {
    uint64_t x_frac = x.FRAC + ( FRAC_SEC_TIME * x.tv_sec );
    uint64_t y_frac = y.FRAC + ( FRAC_SEC_TIME * y.tv_sec );

    x_frac -= y_frac;
    res.tv_sec = x_frac / FRAC_SEC_TIME;
    res.FRAC = x_frac % FRAC_SEC_TIME;

    return x.tv_sec < y.tv_sec;
}

int sum_TIME( TIME &res, TIME &x, TIME &y ) {
    uint64_t x_frac = x.FRAC + ( FRAC_SEC_TIME * x.tv_sec );
    uint64_t y_frac = y.FRAC + ( FRAC_SEC_TIME * y.tv_sec );

    x_frac += y_frac;

    res.tv_sec = x_frac / FRAC_SEC_TIME;
    res.FRAC = x_frac % FRAC_SEC_TIME;

    return 0;
}

int avg_TIME( TIME &avg, TIME &tot, int samples ) {
    uint64_t tmp = tot.FRAC + ( FRAC_SEC_TIME * tot.tv_sec );
    tmp /= samples;

    avg.tv_sec = tmp / FRAC_SEC_TIME;
    avg.FRAC = tmp % FRAC_SEC_TIME;

    return 0;
}

int PrintTime( TIME &t ) {
    AdjustTime( t );
    return printf( TIME_PRINT, t.tv_sec, t.FRAC );
}

void AdjustTime( TIME &t ) {
    uint64_t x_frac = t.FRAC + ( FRAC_SEC_TIME * t.tv_sec );

    t.tv_sec = x_frac / FRAC_SEC_TIME;
    t.FRAC = x_frac % FRAC_SEC_TIME;
}

timespec convertTimeToTimespec( TIME &t ) {
    timespec _t;
    _t.tv_sec = t.tv_sec;
    #if NANO_TIME
        _t.tv_nsec = t.tv_nsec;
    #else
        _t.tv_nsec = t.tv_usec * 1000;
    #endif
    return _t;

}

}
