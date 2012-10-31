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
#ifndef INTERVAL_H
#define INTERVAL_H

#include <fstream>

#include "common.h"

using namespace std;

namespace util {

/**
    An Interval is an abstract class which represents
    a sequential set of positive integers
*/
class Interval {
    public:
        Interval(uint a, uint b);

        uint getStart() const { return start; }
        uint getEnd() const { return end; }

        virtual uint length() = 0;

        virtual ~Interval();
    protected:
        uint start, end;
    private:
};

/**
    An open interval is (a, b) = {x | a < x < b}
*/
class OpenInterval : public Interval {
    public:
        OpenInterval(uint a, uint b) : Interval(a, b) { len = b - a - 1; }
        inline uint length() { return len; }
    protected:
        uint len;
};

/**
    A closed interval is [a, b] = {x | a <= x <= b }
*/
class ClosedInterval : public Interval {
    public:
        ClosedInterval( uint a, uint b) : Interval( a, b) { len = end - start + 1;}
        inline uint length() { return len; }
    protected:
        uint len;
};

/**
    A Left-Closed, Right-Open Interval is [a, b) = { x | a <= x < b }
*/
class LCROInterval : public Interval {
    public:
        LCROInterval( uint a, uint b) : Interval( a, b) { len = end - start; }
        inline uint length() { return len; }

        friend ostream& operator<<(ostream& os, const LCROInterval &interval);
    protected:
        uint len;
};

/**
    A Left-Open, Right-Closed Interval is (a, b] = { x | a < x <= b }
*/
class LORCInterval : public Interval {
    public:
        LORCInterval( uint a, uint b) : Interval( a, b) { len = end - start; }
        inline uint length() { return len; }

        friend ostream& operator<<(ostream& os, const LORCInterval &interval);
    protected:
        uint len;
};

ostream& operator<<(ostream& os, const LCROInterval &interval);
ostream& operator<<(ostream& os, const LORCInterval &interval);

}

#endif // INTERVAL_H
