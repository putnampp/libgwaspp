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
#ifndef TIMING_H_INCLUDED
#define TIMING_H_INCLUDED

#include <sys/time.h>
#include <stdint.h>
#include <cstdio>

#include "common.h"

namespace util {

#ifndef NANO_TIME
#define NANO_TIME 0
#endif

#if NANO_TIME

#define FRAC_SEC_TIME 1000000000
#define FRAC tv_nsec
#define TIME_PRINT "%lu.%09lu s"
#define TIME_ABRV "ns"

#define GetTime(x) clock_gettime(CLOCK_MONOTONIC, &x)
typedef timespec TIME;

#else

#define FRAC_SEC_TIME 1000000
#define FRAC tv_usec
#define TIME_PRINT "%lu.%06lu s"
#define TIME_ABRV "us"

#define GetTime(x) gettimeofday(&x, NULL)
typedef timeval TIME;

#endif

#ifndef INIT_LAPSE_TIME
#if DEBUG_LEVEL
#define INIT_LAPSE_TIME TIME t1, t2, res
#else
#define INIT_LAPSE_TIME
#endif
#endif

#ifndef RECORD_START
#if DEBUG_LEVEL
#define RECORD_START GetTime(t1)
#else
#define RECORD_START
#endif
#endif

#ifndef RECORD_STOP
#if DEBUG_LEVEL
#define RECORD_STOP GetTime(t2)
#else
#define RECORD_STOP
#endif
#endif

#ifndef COMPUTE_LAPSE
#if DEBUG_LEVEL
#define COMPUTE_LAPSE diff_TIME(res, t2, t1);
#else
#define COMPUTE_LAPSE
#endif
#endif

#ifndef PRINT_LAPSE
#if DEBUG_LEVEL
#define PRINT_LAPSE(out, x) COMPUTE_LAPSE; \
                        out << x; \
                        PrintTime(res); \
                        out << "\n";
#else
#define PRINT_LAPSE(out, x)
#endif
#endif


//#define PrintTime(x) printf(TIME_PRINT, x.tv_sec, x.FRAC)

int diff_TIME( TIME &res, TIME &x, TIME &y );
int sum_TIME( TIME &res, TIME &x, TIME &y );
int avg_TIME( TIME &avg, TIME &tot, int samples );

void AdjustTime(TIME &t);
int PrintTime(TIME &t);
timespec convertTimeToTimespec( TIME &t);

}

#endif // TIMING_H_INCLUDED
