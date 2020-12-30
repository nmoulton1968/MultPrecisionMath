/******************************************************************************
Calculate the Natual logarithm e using an infinite series method.
Copyright (C) 1997-2020 Norm Moulton

This is an example program that exercises the MPIM multi-precision integer
class.  It calculates digits of e using an iterative sequence of calculations.
The solution oscillates above and below the correct value and becomes more
precise the longer the program runs.


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

#include "mpim.h"
#include <string.h>

int main()
{
    cout << "Calculating . . .\n";
    cout.flush();

    MPI mCurr;
    MPI mLast;
    MPI mNFact;
    MPI mOffset;
    MPI mTerm;
    int n;
    int j;
    char szCurr[5000]; // decimal representation current computed value
    char szLast[5000]; // decimal representation last computed value
    char *p1, *p2;
    enum
    {
        OFFSET = 1200
    }; // number of decimal digits in offset, (adjust this as needed)

    memset(szCurr, 0, sizeof(szCurr));
    memset(szLast, 0, sizeof(szCurr));
    n = 1;

    mOffset = 1;
    for(j = 0; j < OFFSET; ++j)
        mOffset *= 10;

    mNFact = 1;
    mCurr = mOffset;

    do
    {
        // update last
        mLast = mCurr;
        strcpy(szLast, szCurr);

        // calc next
        mNFact *= n;
        ++n;
        mTerm = mOffset / mNFact;
        mCurr += mTerm;

        mCurr.String(szCurr);

        if(!(n % 20) || mLast == mCurr)
        {
            // how many chars are the same as last time?
            p1 = szCurr;
            p2 = szLast;
            j = 0;
            while(*p1++ == *p2++ && *p1)
                ++j;

            // output current results
            cout << "Terms=" << n << ", Digits=" << j << ", " << mCurr << "\n";
            cout.flush();
        }
    }
    while(mLast != mCurr);

    cout << "Reached limit of calculation capability.\n";
    return 0;
}
