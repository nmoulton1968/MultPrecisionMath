/******************************************************************************
MPIM - Multi Precision Integer Math
Copyright (C) 1997-2020 Norm Moulton

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

#ifndef MPIM_H
#define MPIM_H

// Constants.
#define BASE 10              // Output display base.
#define MPI_BUFF 5000        // Number of chars in decimal representation.
#define INT32 long           // Internal digit storage.
#define INT64 long long      // Intermediate results of digit math.
#define MOD_VALUE 0x40000000 // Largest value stored in a single digit.
#define SHIFT_VALUE 30       // How many bits in MOD_VALUE.
//#define MOD_VALUE 0x100    // Largest value stored in a single digit.
//#define SHIFT_VALUE 8      // How many bits in MOD_VALUE.

// Number of digits in internal representation.
#define MAX_ARRAY 500 // 15000 bits.

class MPI
{
public:  // Data.
    INT32 mArray[MAX_ARRAY];
    bool mIsOverflow;

public: // Functions.

    // Constructors/ Assignments.
    MPI();
    void Zero();  // Set to zero.
    MPI(char*);   // Construct from a decimal string.
    MPI(int);     // Construct from an integer.

    MPI operator=(const MPI&);
    MPI operator=(int);
    MPI operator=(char*);

    // Arithmetic, eg.  x = y + z.
    MPI operator+(const MPI&) const;
    MPI operator+(int) const;
    MPI operator-(const MPI&) const;
    MPI operator-(int) const;
    MPI operator*(const MPI&) const;
    MPI operator*(int) const;
    MPI operator/(const MPI&) const;
    MPI operator/(int) const;
    MPI operator%(const MPI&) const;
    MPI operator%(int) const;
    MPI operator^(const MPI&) const;
    MPI operator^(int) const;

    MPI operator++();
    MPI operator++(int);
    MPI operator--();
    MPI operator--(int);

    // Shortcut Forms, eg. x += y.
    MPI operator+=(const MPI&);
    MPI operator+=(int);
    MPI operator-=(const MPI&);
    MPI operator-=(int);
    MPI operator*=(const MPI&);
    MPI operator*=(int);
    MPI operator/=(const MPI&);
    MPI operator/=(int);
    MPI operator%=(const MPI&);
    MPI operator%=(int);
    MPI operator^=(const MPI&);
    MPI operator^=(int);

    // Multiplication: Simple method.
    MPI MultSmpl(const MPI& y) const;

    // Multiplication: Divide and Conquer.
    MPI MultDC(const MPI&) const;

    // Multiplication: A La Russe.
    MPI MultALR(const MPI&) const;

    // Special Divide: Return quotient and remainder.
    MPI Divide(const MPI&, MPI&) const;

    // Modulus Functions.
    MPI ModMult(const MPI&, const MPI&) const;
    MPI ModPow(const MPI&, const MPI&) const;

    // Comparison/ Logical, eg. if(x < y).
    bool operator<(const MPI&) const;
    bool operator>(const MPI&) const;
    bool operator<=(const MPI&) const;
    bool operator>=(const MPI&) const;
    bool operator==(const MPI&) const;
    bool operator!=(const MPI&) const;

    // Conversions and I/O.
    int Integer() const;          // Convert to integer.
    char* String(char []) const;  // Convert to decimal char buffer.
    friend ostream& operator<<(ostream&, MPI&);
    friend istream& operator>>(istream&, MPI&);

    // Helper Functions and Diagnostics
    bool IsValid() const;                 // Check validity of MPI.
    void Display() const;                 // Show internal representation.
    inline void ShiftLeft(const int);     // Shift internal digits n positions.
    inline void ShiftRight(const int);    // Shift internal digits n positions.
    inline void Mult2();                  // Quick multiply (using shift) by 2.
    inline void Div2();                   // Quick divide (using shift) by 2.

    INT32 MSDigit() const;                // Most significant digit.
    int Size() const;                     // Count number of significant digits.
    inline int Largest(const MPI&) const; // Return Size() of largest.

};
#endif

