/******************************************************************************
MPIM - Multi Precision Integer Math
Copyright (C) 1997-2020 Norm Moulton

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
#include <cstring>

/*****************************************************************************/
// CONSTRUCTORS and CONVERSIONS
/*****************************************************************************/

// Construct and initialize to zero.
MPI::MPI()
{
    Zero();
}

// Initialize value to zero.
void MPI::Zero()
{
    int i;
    for(i=0; i<MAX_ARRAY; ++i)
    {
        mArray[i] = 0;
    }

    mIsOverflow = false;
}

// Construct from decimal string representation.
MPI::MPI(char* psz)
{
    // Validate input.
    if(psz == 0)
    {
        return;
    }

    Zero();

    // For each digit . . .
    while(*psz != 0)
    {
        // Multiply by ten.
        *this = *this * 10;

        // Add digit value.
        int digit = *psz - '0';
        *this = *this + digit;

        // Move to next digit.
        ++psz;
    }
}

// Construct from integer.
MPI::MPI(int n)
{
    Zero();
    mArray[0] = n & (MOD_VALUE-1);
}

/*****************************************************************************/
// ASSIGNMENTS
/*****************************************************************************/

MPI MPI::operator=(const MPI& m)
{
    // Protect against assigning object to itself.
    if(&m == this) return *this;

    // Duplicate the array.
    for(int i=0; i<MAX_ARRAY; ++i)
    {
        mArray[i] = m.mArray[i];
    }

    mIsOverflow = m.mIsOverflow;

    return *this;
}

MPI MPI::operator=(int n)
{
    *this = MPI(n);
    return *this;
}

MPI MPI::operator=(char* psz)
{
    *this = MPI(psz);
    return *this;
}

/*****************************************************************************/
// ADDITION AND SUBTRACTION
/*****************************************************************************/

// Add MPI + MPI.
// Algorithm based on Menezes, 14.7, p. 594.
MPI MPI::operator+(const MPI& m) const
{
    MPI w;    // result
    int carry = 0;

    for(int i=0; i<MAX_ARRAY; ++i)
    {
        // Add current digits.
        w.mArray[i] = mArray[i] + m.mArray[i] + carry;

        // Test for digit overflow.
        if(w.mArray[i] & MOD_VALUE)
        {
            w.mArray[i] -= MOD_VALUE;
            carry = 1;
        }
        else
        {
            carry = 0;
        }
    }

    // If still carry, whole result is overflow.
    w.mIsOverflow |= (carry != 0);

    return w;
}

// Add MPI + int.
// Algorithm based on Menezes, 14.7, p. 594.
MPI MPI::operator+(int n) const
{
    MPI w;

    n = n & (MOD_VALUE-1);

    // Add argument digit.
    w.mArray[0] = mArray[0] + n;

    // Digit overflow?
    int carry = 0;
    if(w.mArray[0] & MOD_VALUE)
    {
        w.mArray[0] -= MOD_VALUE;
        carry = 1;
    }

    // Copy all other digits, propagating possible carry.
    for(int i=1; i<MAX_ARRAY; ++i)
    {
        // Copy digit, add carry.
        w.mArray[i] = mArray[i] + carry;

        // Test for digit overflow.
        if(w.mArray[i] & MOD_VALUE)
        {
            w.mArray[i] -= MOD_VALUE;
            carry = 1;
        }
        else
        {
            carry = 0;
        }
    }

    // If carry still exists, whole result is an overflow.
    w.mIsOverflow |= (carry != 0);

    return w;
}

// Subtract MPI - MPI.
// Algorithm based on Menezes, 14.9, p. 595.
MPI MPI::operator-(const MPI& m) const
{
    MPI w;
    int carry = 0;

    for(int i=0; i<MAX_ARRAY; ++i)
    {
        // Subtract current digits.
        w.mArray[i] = mArray[i] - m.mArray[i] - carry;

        // Test for digit overflow.
        if(w.mArray[i] & MOD_VALUE)
        {
            w.mArray[i] += MOD_VALUE;
            carry = 1;
        }
        else
        {
            carry = 0;
        }
    }

    // If carry still exists, whole result is an overflow.
    w.mIsOverflow |= (carry != 0);

    return w;
}

// Subtract MPI - int.
// Algorithm based on Menezes, 14.9, p. 595.
MPI MPI::operator-(int n) const
{
    MPI w;

    n = n & (MOD_VALUE-1);

    // Subtract argument digit.
    w.mArray[0] = mArray[0] - n;

    // Digit overflow?
    int carry = 0;
    if(w.mArray[0] & MOD_VALUE)
    {
        w.mArray[0] += MOD_VALUE;
        carry = 1;
    }

    for(int i=1; i<MAX_ARRAY; ++i)
    {
        // Subtract current digits.
        w.mArray[i] = mArray[i] - carry;

        // Test for digit overflow.
        if(w.mArray[i] & MOD_VALUE)
        {
            w.mArray[i] += MOD_VALUE;
            carry = 1;
        }
        else
        {
            carry = 0;
        }
    }

    // If carry still exists, whole result is an overflow.
    w.mIsOverflow |= (carry != 0);

    return w;
}

/*****************************************************************************/
// MULTIPLICATION
/*****************************************************************************/

#define BREAK_EVEN 150

// Recursize multiply controller.
MPI MPI::operator*(const MPI& y) const
{
    if(Size() < BREAK_EVEN || y.Size() < BREAK_EVEN)
    {
        return MultSmpl(y);
    }
    else
    {
        return MultDC(y);
    }
}

// Multiply MPI * MPI.
// Algorithm based on Menezes, 14.12, p. 595.
MPI MPI::MultSmpl(const MPI& y) const
{
    MPI w;  // result

    INT64 uv; // double the precision of single digit
    INT64 carry;

    int n = Size();   // Get # digits in x.
    int t = y.Size(); // Get # digits in y.

    for(int i=0; i<t; ++i)
    {
        carry = 0;

        int j;
        for(j=0; j<n; ++j)
        {
            if(i+j > MAX_ARRAY-1)
                continue;

            uv = (INT64)w.mArray[i+j] + (INT64)mArray[j] * (INT64)y.mArray[i] + carry;
            w.mArray[i+j] = (INT32)(uv & (MOD_VALUE-1));  // LOWORD.
            carry = ((uv >> SHIFT_VALUE) & (MOD_VALUE-1));  // HIWORD.
        }

        if(i+j < MAX_ARRAY)
        {
            w.mArray[i+j] += (INT32)carry;
        }
        else if(carry > 0)
        {
            w.mIsOverflow = true;
        }
    }

    return w;
}

// Multiply MPI * int.
// Algorithm based on Menezes, 14.12, p. 595.
MPI MPI::operator*(int y) const
{
    MPI w;  // result
    int n;  // # digits in x

    INT64 uv; // double the precision of single digit
    INT64 carry;

    n = Size()+1;   // Get # digits in x + 1.
    y = y & (MOD_VALUE-1);
    carry = 0;

    for(int j=0; j<n; ++j)
    {
        uv = (INT64)mArray[j] * (INT64)y + carry;
        w.mArray[j] = (INT32)(uv & (MOD_VALUE-1));   // LOWORD
        carry = ((uv >> SHIFT_VALUE) & (MOD_VALUE-1)); // HIWORD
    }

    if(carry > 0)
    {
        w.mIsOverflow = true;
    }

    return w;
}

// Multiply MPI * MPI, A la Russe.
// Algorithm based on Brassard, p. 4.
MPI MPI::MultALR(const MPI& m) const
{
    MPI w; // result
    MPI x; // multiplicand
    MPI y; // multiplier

    // Optimize order; better when multiplier is smaller.
    if(x < y)
    {
        x = m;
        y = *this;
    }
    else
    {
        x = *this;
        y = m;
    }

    // Calculate how many loops to perform.
    int j = x.Size() * SHIFT_VALUE;

    for(int i=0; i<j; ++i)
    {
        // Sum if odd.
        if(y.mArray[0] & 1)
        {
            w += x;
        }

        x.Mult2();
        y.Div2();
    }

    return w;
}

// Multiply MPI * MPI, Divide and Conquer.
// Algorithm based on Brassard, p. 219-223.
MPI MPI::MultDC(const MPI& m) const
{
    MPI d;           // result product
    MPI w, x, y, z;  // split pieces of the arguments
    MPI p, q, r;     // intermediate products
    int f, h;        // full and half size of arguements

    f = Largest(m)-1;
    h = f/2;

    // Assign w.
    for(int i=h; i<f; ++i)
    {
        w.mArray[i-h] = mArray[i];
    }

    // Assign x.
    for(int i=0; i<h; ++i)
    {
        x.mArray[i] = mArray[i];
    }

    // Assign y.
    for(int i=h; i<f; ++i)
    {
        y.mArray[i-h] = m.mArray[i];
    }

    // Assign z.
    for(int i=0; i<h; ++i)
    {
        z.mArray[i] = m.mArray[i];
    }

    // Calculate three products.
    p = w * y;
    q = x * z;
    r = (w + x) * (y + z);

    d = p;
    d.ShiftLeft(h);

    d += r - p - q;
    d.ShiftLeft(h);

    d += q;

    return d;
}

/*****************************************************************************/
// DIVISION AND MODULUS
/*****************************************************************************/

// Divide MPI / MPI, Classical.
// Algorithm based on Knuth, D, p. 257.
MPI MPI::operator/(const MPI& m) const
{
    MPI u;    // dividend
    MPI v;    // divisor
    MPI q;    // quotient
    MPI s;    // trial subtract amount
    INT32 qh; // trial quotient
    INT32 d;  // normalization factor

    // Work with local copies
    u = *this;
    v = m;

#if 0
    if(u.Size() == 0)
        return MPI(0);

    if(v.Size() == 0)
    {
        q.mIsOverflow = true;
        return q;
    }
#endif

    // Normalize.
    d = MOD_VALUE / (v.MSDigit() + 1);
    u *= (int)d;
    v *= (int)d;

    // Get # digits.
    int t = v.Size();
    int n = u.Size();

    // Main calculation loop.
    for(int j=n; j>=t; --j)
    {
        // Calculate trial quotient.
        qh = (((INT64)u.mArray[j] * MOD_VALUE) + u.mArray[j-1]) / v.mArray[t-1];

        // Adjust quotient if too large.
        if(qh & MOD_VALUE)
            qh = MOD_VALUE -1;

        // Multiply and subtract.
        s = v;
        s.ShiftLeft(j-t);
        u -= s * qh;

        // Adjust result if it went negative.
        while(u.mArray[MAX_ARRAY-1] == MOD_VALUE-1)
        {
            --qh;   // Adjust quotient digit.
            u += s; // Add back.
            u.mIsOverflow = false;
        }

        // Set the quotient digit we just found.
        q.mArray[j-t] = qh;
    }

    return q;
}

// Divide MPI / int.
// Algorithm based on Knuth, D, p. 257.
MPI MPI::operator/(int y) const
{
    MPI u;    // dividend
    INT32 v;  // divisor
    MPI q;    // quotient
    MPI s;    // trial subtract amount
    INT32 qh; // trial quotient
    INT32 d;  // normalization factor

    // Work with local copies.
    u = *this;
    v = y & (MOD_VALUE-1);

    // Normalize.
    d = MOD_VALUE / (v + 1);
    u *= (int)d;
    v *= (int)d;

    // Get # digits.
    int n = u.Size(); // # digits in u.

    // Main calculation loop.
    for(int j=n; j>=1; --j)
    {
        // Calculate trial quotient.
        qh = (((INT64)u.mArray[j] * MOD_VALUE) + u.mArray[j-1]) / v;

        // Adjust quotient if too large.
        if(qh & MOD_VALUE)
            qh = MOD_VALUE -1;

        // Multiply and subtract.
        s = v;
        s.ShiftLeft(j-1);
        u -= s * qh;

        // Adjust result if it went negative.
        while(u.mArray[MAX_ARRAY-1] == MOD_VALUE-1)
        {
            --qh;   // Adjust quotient digit.
            u += s; // Add back.
            u.mIsOverflow = false;
        }

        // Set the quotient digit we just found.
        q.mArray[j-1] = qh;
    }

    return q;
}

// Modulus, MPI % MPI.
// Algorithm based on Knuth, D, p. 257.
MPI MPI::operator%(const MPI& m) const
{
    MPI u;    // dividend
    MPI v;    // divisor
    MPI s;    // trial subtract amount
    INT32 qh; // trial quotient
    INT32 d;  // normalization factor

    // Work with local copies.
    u = *this;
    v = m;

    // Normalize.
    d = MOD_VALUE / (v.MSDigit() + 1);
    u *= (int)d;
    v *= (int)d;

    // Get # digits.
    int t = v.Size(); // # digits in u.
    int n = u.Size(); // # digits in v.

    // Main calculation loop.
    for(int j=n; j>=t; --j)
    {
        // Calculate trial quotient, get LOWORD.
        qh = (((INT64)u.mArray[j] * MOD_VALUE) + u.mArray[j-1]) / v.mArray[t-1];

        // Adjust quotient if too large.
        if(qh & MOD_VALUE)
            qh = MOD_VALUE -1;

        // Multiply and subtract.
        s = v;
        s.ShiftLeft(j-t);
        u -= s * qh;

        // Adjust result if it went negative.
        while(u.mArray[MAX_ARRAY-1] == MOD_VALUE-1)
        {
            u += s; // Add back.
        }
    }

    // Un-normalize.
    u /= d;

    // Fix flag if it ever went negative.
    u.mIsOverflow = false;

    return u;  // The remainder.
}

// Modulus MPI % int.
MPI MPI::operator%(int n) const
{
    return *this % MPI(n);
}

// Division, Quotient and Remainder, MPI / MPI.
// Algorithm based on Knuth, D, p. 257.
MPI MPI::Divide(const MPI& v1, MPI& u) const
{
    MPI q;    // quotient
    MPI v;    // divisor
    MPI s;    // trial subtract amount
    INT32 qh; // trial quotient
    INT32 d;  // normalization factor

    // Work with local copies.
    u = *this;
    v = v1;

    // Normalize.
    d = MOD_VALUE / (v.MSDigit() + 1);
    u *= (int)d;
    v *= (int)d;

    // Get # digits.
    int t = v.Size(); // # digits in u.
    int n = u.Size(); // # digits in v.

    // Main calculation loop.
    for(int j=n; j>=t; --j)
    {
        // Calculate trial quotient, get LOWORD.
        qh = (((INT64)u.mArray[j] * MOD_VALUE) + u.mArray[j-1]) / v.mArray[t-1];

        // Adjust quotient if too large.
        if(qh & MOD_VALUE)
            qh = MOD_VALUE -1;

        // Multiply and subtract.
        s = v;
        s.ShiftLeft(j-t);
        u -= s * qh;

        // Adjust result if it went negative.
        while(u.mArray[MAX_ARRAY-1] == MOD_VALUE-1)
            u += s; // Add back.

        // Set the quotient digit we just found.
        q.mArray[j-t] = qh;
    }

    // Un-normalize.
    u /= d;

    // Fix flag if it ever went negative.
    u.mIsOverflow = false;

    // Remainder is returned in u.
    return q;
}

/*****************************************************************************/
// EXPONENTS AND MODULAR ARITHMETIC
/*****************************************************************************/

// Exponential MPI ^ MPI.
// Algorithm based on CLR, p. 829.
MPI MPI::operator^(const MPI& y) const
{
    MPI w;  // return value
    MPI s;  // shifted exponent
    int n;  // digits in exponent
    int k;  // bits in exponent

    // shifted exponent
    s = y;

    // How many times to loop.
    n = y.Size();
    k = n * SHIFT_VALUE;

    // x ^ 0 = 1.
    if(!n) return MPI(1);

    // Main loop.
    w = 1;
    for(int i=0; i<k; ++i)
    {
        // Square.
        w *= w;

        // Multiply.
        if(s.mArray[n-1] & MOD_VALUE/2)
            w *= *this;

        // Shift.
        s.Mult2();
    }

    return w;
}

// Exponential MPI ^ int.
MPI MPI::operator^(int n) const
{
    return *this ^ MPI(n);
}

// Modular Multiplication.
// Algorithm based on Menezes, 14.28 p. 600.
MPI MPI::ModMult(const MPI& y1, const MPI& m) const
{
    MPI w;
    MPI x;
    MPI y;

    x = *this;
    y = y1;

    x %= m;
    y %= m;

    w = x * y;
    w %= m;

    return w;
}

// Modular Exponetial, Repeated Squaring.
// Algorithm based on CLR, p. 829.
MPI MPI::ModPow(const MPI& y, const MPI& m) const
{
    MPI w;  // return value
    MPI s;  // shifted exponent
    int n;  // digits in exponent
    int k;  // bits in exponent

    // Shifted exponent.
    s = y;

    // How many times to loop.
    n = y.Size();
    k = n * SHIFT_VALUE;

    // x ^ 0 = 1.
    if(!n) return MPI(1);

    // Main loop.
    w = 1;
    for(int i=0; i<k; ++i)
    {
        // Square.
        w *= w;
        w %= m;

        // Multiply.
        if(s.mArray[n-1] & MOD_VALUE/2)
        {
            w *= *this;
            w %= m;
        }

        // Shift.
        s.Mult2();
    }

    return w;
}

/*****************************************************************************/
// ARITHMETIC SHORTCUT FORMS
/*****************************************************************************/

MPI MPI::operator+=(const MPI& m)
{
    return *this = *this + m;
}

MPI MPI::operator-=(const MPI& m)
{
    return *this = *this - m;
}

MPI MPI::operator*=(const MPI& m)
{
    return *this = *this * m;
}

MPI MPI::operator/=(const MPI& m)
{
    return *this = *this / m;
}

MPI MPI::operator%=(const MPI& m)
{
    return *this = *this % m;
}

MPI MPI::operator^=(const MPI& m)
{
    return *this = *this ^ m;
}

MPI MPI::operator+=(int n)
{
    return *this = *this + n;
}

MPI MPI::operator-=(int n)
{
    return *this = *this - n;
}

MPI MPI::operator*=(int n)
{
    return *this = *this * n;
}

MPI MPI::operator/=(int n)
{
    return *this = *this / n;
}

MPI MPI::operator%=(int n)
{
    return *this = *this % n;
}

MPI MPI::operator^=(int n)
{
    return *this = *this ^ n;
}

MPI MPI::operator++()
{
    *this = *this + 1;
    return *this;
}
MPI MPI::operator++(int)
{
    MPI m = *this;
    *this = *this + 1;
    return m;
}
MPI MPI::operator--()
{
    *this = *this - 1;
    return *this;
}
MPI MPI::operator--(int)
{
    MPI m = *this;
    *this = *this - 1;
    return m;
}

/*****************************************************************************/
// LOGICAL OPERATORS/ COMPARISON/ SHIFTING
/*****************************************************************************/

bool MPI::operator<(const MPI& m) const
{
    MPI x;

    x = *this - m;
    return x.mIsOverflow;
}

bool MPI::operator<=(const MPI& m) const
{
    MPI x;

    x = m - *this;
    return !x.mIsOverflow;
}

bool MPI::operator>(const MPI& m) const
{
    MPI x;

    x = m - *this;
    return x.mIsOverflow;
}

bool MPI::operator>=(const MPI& m) const
{
    MPI x;

    x = *this - m;
    return !x.mIsOverflow;
}

bool MPI::operator==(const MPI& m) const
{
    for(int i=0; i<MAX_ARRAY; ++i)
    {
        if(m.mArray[i] != mArray[i])
        {
            return false;
        }
    }

    return true;
}

bool MPI::operator!=(const MPI& m) const
{
    for(int i=0; i<MAX_ARRAY; ++i)
    {
        if(m.mArray[i] != mArray[i])
        {
            return true;
        }
    }

    return false;
}

/*****************************************************************************/
// CONVERSIONS AND I/O
/*****************************************************************************/

// Convert to integer.
int MPI::Integer() const
{
    // If the value is too large, flag an error.
    if(Size() > 1)
    {
        return -1;
    }

    return (int)mArray[0];
}

// Convert to a character string representation in decimal.
char* MPI::String(char sz[]) const
{
    MPI q;       // quotient
    MPI r;       // remaider
    MPI zero;
    MPI base;

    // Don't print invalid number.
    if(Size() >= MAX_ARRAY || mIsOverflow)
    {
        strcpy(sz, "ERROR");
        return sz;
    }

    base = BASE;
    int i = 0;
    q = *this;

    do
    {
        q = q.Divide(base, r);
        sz[i++] = r.mArray[0] + '0';
    }
    while(q != zero);

    sz[i] = '\0';
    --i;

    // Reverse the string into human readable order.
    char temp;
    for(int j=0; j<=(i/2); ++j)
    {
        temp = sz[j];
        sz[j] = sz[i-j];
        sz[i-j] = temp;
    }

    return sz;
}

// Stream output.
ostream& operator<<(ostream& os, MPI& m)
{
    char sz[MPI_BUFF];

    m.String(sz);
    os << sz;
    return os;
}

// Stream input.
istream& operator>>(istream& is, MPI& m)
{
    char sz[MPI_BUFF];

    is >> sz;
    m = MPI(sz);
    return is;
}

/*****************************************************************************/
// HELPER FUNCTIONS AND DIAGNOSTICS
/*****************************************************************************/

// Check for digit overflow and overflow flag.
bool MPI::IsValid() const
{
    if(mIsOverflow)
    {
        return false;
    }

    for(int i=0; i<MAX_ARRAY; ++i)
    {
        if(mArray[i] > MOD_VALUE || mArray[i] < 0)
        {
            return false;
        }
    }

    return true;
}

// Diagnostic display of the internal representation.
void MPI::Display() const
{
    cout << "MPI [";
    for(int i=0; i<MAX_ARRAY; ++i)
    {
        cout << mArray[i] << ", ";
    }

    cout << "]\n";
}

// Most significant digit.
INT32 MPI::MSDigit() const
{
    for(int i=MAX_ARRAY-1; i>=0; --i)
    {
        if(mArray[i] != 0)
        {
            return mArray[i];
        }
    }

    return 0;
}

// Count the number of significant digits.
int MPI::Size() const
{
    for(int i=MAX_ARRAY-1; i>=0; --i)
    {
        if(mArray[i] != 0)
        {
            return i+1;
        }
    }

    return 0;
}

// Return the largest size of two arguments.
inline int MPI::Largest(const MPI& m) const
{
    int a, b, s;

    a = Size();   // Size of left argument.
    b = m.Size(); // Size of right argument.

    // Find larger argument.
    if(a > b)
    {
        s = a;
    }
    else
    {
        s = b;
    }

    return s + 1;
}

// Shift by n digit positions.
void MPI::ShiftLeft(const int n)
{
    if(n >= MAX_ARRAY-1)
    {
        mIsOverflow = true;
        return;
    }

    for(int i=MAX_ARRAY-1; i>n-1; --i)
    {
        mArray[i] = mArray[i-n];
    }

    for(int i=n-1; i>=0; --i)
    {
        mArray[i] = 0;
    }
}

// Shift by n digit positions.
void MPI::ShiftRight(const int n)
{
    if(n >= MAX_ARRAY-1)
    {
        mIsOverflow = true;
        return;
    }

    for(int i=n; i<MAX_ARRAY; ++i)
    {
        mArray[i-n] = mArray[i];
    }

    for(int i=MAX_ARRAY-n; i<MAX_ARRAY; ++i)
    {
        mArray[i] = 0;
    }
}

// Simple multiply using bit shift.
void MPI::Mult2()
{
    INT32 carry = 0;

    for(int i=0; i<MAX_ARRAY; ++i)
    {
        mArray[i] <<= 1;
        if(carry) mArray[i] |= 1;

        carry = mArray[i] & MOD_VALUE;
        mArray[i] &= MOD_VALUE-1;
    }
}

// Simple divide using bit shift.
void MPI::Div2()
{
    INT32 carry = 0;

    for(int i=MAX_ARRAY-1; i>=0; --i)
    {
        if(carry) mArray[i] |= MOD_VALUE;

        carry = mArray[i] & 1;
        mArray[i] >>= 1;
    }
}
