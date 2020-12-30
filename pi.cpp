/******************************************************************************
Calculate Pi using infinite series arctangent method
Copyright (C) 1997-2020 Norm Moulton

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

#include <cstring>
#include "mpim.h"


class ArcTan
{
private:
    MPI mX;
    MPI mValue;
    MPI mOffset;
    int mExponent;
    enum { OFFSET = 1000 };
    const MPI ZERO;

public:
    // Ctor from int
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
        MPI iTerm;
        iTerm = (mX ^ mExponent);
        iTerm *= mOffset;
        iTerm /= mExponent;

        for(int j=1; j<mExponent; ++j)
        {
            iTerm /= 1000;
        }

        if(iTerm == ZERO)
            return mValue;

        mValue += iTerm;
        mExponent += 2;

        iTerm = (mX ^ mExponent);
        iTerm *= mOffset;
        iTerm /= mExponent;

        for(int j=1; j< mExponent; ++j)
        {
            iTerm /= 1000;
        }

        if(iTerm == ZERO) return mValue;
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
    ArcTan arcTanOneHalf(500);   // arctan(1/2)
    ArcTan arcTanOneFifth(200);  // arctan(1/5)
    ArcTan arcTanOneEighth(125); // arctan(1/8)
    char sCurr[MPI_BUFF] = {0};
    char sLast[MPI_BUFF] = {0};
    char *p1, *p2;

    // Digits found history.
    int digitsFound1 = 0;
    int digitsFound0 = 0;

    const int REPORT_ITERATIONS = 10;

    bool isDone = false;
    int i = 0;
    do
    {
        ++i;
        mLast = mCurr;
        mCurr = (arcTanOneHalf.Next() + arcTanOneFifth.Next() + arcTanOneEighth.Next())* 4;

        // Periodically output the current estimate.
        if((i%REPORT_ITERATIONS) == 0)
        {
            // Have we exceeded the resolution of the registers?
            if(mLast == mCurr) isDone = true;

            strncpy(sLast, sCurr, sizeof(sLast));
            mCurr.String(sCurr);

            // Count how many chars match in the current and previous estimates.
            // We assume digits that have stablilized are correct digits, which
            // is largely true due to the convergent nature of the algorithm.
            // This might overcount by one or two digits, so it is a rough metric.
            p1 = sCurr;
            p2 = sLast;

            // Update digit found history.
            digitsFound1 = digitsFound0;
            digitsFound0 = 0;

            while(*p1++ == *p2++)
            {
                ++digitsFound0;
            }

            // The solution should converge at a roughly constant rate.
            // This rate has been seen to be about 1.1 to 1.4 digits per iteration.
            // When the rate of finding digits increases signifigantly, (more than 2
            // per iteration), it indicates we have probably exceeded the accuracy
            // achievable with the current OFFSET value and must stop. Any further
            // digits after this point will not be accurate.
            int digitRate = digitsFound0 - digitsFound1;
            if((digitRate >= 0) && ((digitRate == 0) || (digitRate < REPORT_ITERATIONS * 2)))
            {
                // Copy to the display buffer only the valid chars found.
                char sDisp[MPI_BUFF] = {0};
                strncpy(sDisp, sCurr, digitsFound0);

                // Output current result.
                cout << "Terms=" << i*2 << ", Digits=" << digitsFound0;
                cout << ", Rate=" << digitRate << "\n";
                cout << sDisp << "\n";
                cout.flush();
            }
            else
            {
                isDone = true;
                cout << "[" << digitsFound1 << ", " << digitsFound0 << "]\n";
            }
        }
    }
    while(!isDone);

    cout << "Reached limit of calculation capability.\n";

    return 0;
}

