/******************************************************************************
Calculate Pi using infinite series arctangent method
Copyright (C) 1997-2016 Norm Moulton

This is an example program that exercises the MPIM multi-precision integer
class.  It calculates digits of Pi using an iterative sequence of arctangent
calculations.  The solution oscillates above and below the correct value, and
becomes more precise the longer the program runs.


This program is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

******************************************************************************/

#include <string.h>
#include "mpim.h"


class ArcTan
{
private:
    MPI mX;
    MPI mValue;
    MPI mOffset;
    int mExponent;
    enum { OFFSET = 1000 };
    const MPI mZero;

public:
    // ctor from int
    ArcTan(int x)
    {
        int j;

        mValue = 0;
        mOffset = 1;
        mExponent = 1;

        mX = x;

        for(j=0; j<OFFSET; ++j)
        {
            mOffset *= 10;
        }
    }

    MPI Curr()
    {
        return mValue;
    }

    MPI Next()
    {
        int j;
        MPI iTerm;
        iTerm = (mX ^ mExponent);
        iTerm *= mOffset;
        iTerm /= mExponent;

        for(j=1; j<mExponent; ++j)
        {
            iTerm /= 1000;
        }

        if(iTerm == mZero)
            return mValue;

        mValue += iTerm;
        mExponent += 2;

        iTerm = (mX ^ mExponent);
        iTerm *= mOffset;
        iTerm /= mExponent;

        for(j=1; j< mExponent; ++j)
        {
            iTerm /= 1000;
        }

        if(iTerm == mZero) return mValue;
        mValue -= iTerm;
        mExponent += 2;

        return mValue;
    }
};

int main(int argc, char* argv[])
{
    cout << "Calculating . . .\n";
    cout.flush();

    MPI mCurr, mLast;
    int i, j;
    ArcTan a500(500); // arctan(1/2)
    ArcTan a200(200); // arctan(1/5)
    ArcTan a125(125); // arctan(1/8)
    char sCurr[5000];
    char sLast[5000];
    char *p1, *p2;

    i = 0;
    sCurr[0] = 0;
    sLast[0] = 0;

    do
    {
        ++i;
        mLast = mCurr;
        mCurr = (a500.Next() + a200.Next() + a125.Next())* 4;

        if(!(i%2))
        {
            strcpy(sLast, sCurr);
            mCurr.String(sCurr);
            // how many chars same as last time?
            p1 = sCurr;
            p2 = sLast;
            j = 0;
            while(*p1++ == *p2++)
            {
                ++j;
            }
            //output current result
            cout << "Terms=" << i*2 << ", Digits=" << j << ", " << mCurr << "\n";
            cout.flush();
        }
    }
    while(mLast != mCurr);

    char c;
    cin >> c;

    return 0;
}

