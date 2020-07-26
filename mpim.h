/******************************************************************************
MPIM - Multi Precision Integer Math
Copyright (C) 1997-2016 Norm Moulton

This class provides a set of functions to perform multi-precision integer
arithmetic operations. It defines the MPI multi-byte numeric type with a
high level of integration into the C++ syntax. Operator overloads and
conversions are provided to support natural expressions combining the
MPI type with native integer operations.

In addition, input and output conversions are provided for ASCII string
conversions, which are useful for initialization, as well as for human
readable display of results.

The class can be configured at compile time to set the size of the native
machine register and the maximum size of the supported MPI type. There is
a trade-off between the supported size and the speed of calculations.


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

#include <iostream>
using namespace std;

// general
#define BOOL int
#define TRUE -1
#define FALSE 0

// MPI attributes
#define MPI_BUFF 20000     // number of chars in decimal representation
#define MAX_ARRAY 5000      // number of digits in internal representation
#define INT32 __int32     // internal digit storage
#define INT64 __int64     // intermediate results of digit math
//#define MOD_VALUE 0x100 // largest value stored in a single digit
//#define SHIFT_VALUE 8   // how many bits in MOD_VALUE
#define MOD_VALUE 0x10000 // largest value stored in a single digit
#define SHIFT_VALUE 16    // how many bits in MOD_VALUE

class MPI
{
public:  // data
    INT32 m_ar[MAX_ARRAY];
    BOOL m_bOverflow;

public: // functions

    // constructors/conversions
    MPI();
    MPI(char* psz); // construct from a decimal string
    MPI(INT32 n);   // construct from an integer
    void Clear();   // set to zero
    MPI operator=(const MPI&);
    MPI operator=(int n);
    MPI operator=(char* psz);


    // normal arithmetic, eg.  x = y + z
    MPI operator+(const MPI& m) const;
    MPI operator+(const int n) const;
    MPI operator-(const MPI& m) const;
    MPI operator-(const int n) const;
    MPI operator*(const MPI& m) const;
    MPI operator*(const int n) const;
    MPI operator/(const MPI& m) const;
    MPI operator/(const int n) const;
    MPI operator%(const MPI& m) const;
    MPI operator%(const int n) const;
    MPI operator^(const MPI& m) const;
    MPI operator^(const int n) const;
    MPI operator++();
    MPI operator++(int);
    MPI operator--();
    MPI operator--(int);

    // shortcut forms, eg. x += y
    MPI operator+=(const MPI& m);
    MPI operator+=(const int n);
    MPI operator-=(const MPI& m);
    MPI operator-=(const int n);
    MPI operator*=(const MPI& m);
    MPI operator*=(const int n);
    MPI operator/=(const MPI& m);
    MPI operator/=(const int n);
    MPI operator%=(const MPI& m);
    MPI operator%=(const int n);
    MPI operator^=(const MPI& m);
    MPI operator^=(const int n);

    // comparison/ logical, eg. if(x < y)
    BOOL operator<(const MPI& m) const;
    BOOL operator>(const MPI& m) const;
    BOOL operator<=(const MPI& m) const;
    BOOL operator>=(const MPI& m) const;
    BOOL operator==(const MPI& m) const;
    BOOL operator!=(const MPI& m) const;

    // Multiplication by Divide and Conquer
    MPI MultDC(MPI& m) const;

    // Multiplication A La Russe
    MPI MultALR(MPI& m) const;

    // modulus versions, eg. x.ModMult(y, m)
    MPI ModMult(const MPI& y, const MPI& m) const;
    MPI ModPow(const MPI& y, const MPI& m) const;

    // c++ stream i/o
    friend ostream& operator<<(ostream& os, MPI& m);
    friend istream& operator>>(istream& is, MPI& m);
    char* String(char a[]) const;  // convert MPI into decimal char buffer

    // diagnostic functions
    BOOL IsValid() const;          // check validity of MPI
    void Display() const;          // show internal representation

    // helper functions
    inline void ShiftLeft(const int n);   // shift internal digits n positions
    inline void ShiftRight(const int n);  // shift internal digits n positions
    inline void Mult2();                  // quick multiply (using shift) by 2
    inline void Div2();                   // quick divide (using shift) by 2

    INT32 MSDigit() const;         // most significant digit
    int Size() const;              // count number of significant digits
    inline int Largest(MPI& m) const;  // return Size() of largest

};

