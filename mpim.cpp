/******************************************************************************
MPIM - Multi Precision Integer Math
Copyright (C) 1997-2016 Norm Moulton

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

/*****************************************************************************/
// CONSTRUCTORS and CONVERSIONS
/*****************************************************************************/

// construct and initialize to zero
MPI::MPI()
{
    Clear();
}

// set value to zero
void MPI::Clear()
{
    int i;
    for(i=0; i<MAX_ARRAY; ++i)
        m_ar[i] = 0;

    m_bOverflow = FALSE;
}

// convert from decimal string representation
MPI::MPI(char* psz)
{
    int nDigit;

    Clear();

    // error
    if(psz == NULL)
        return;

    // for each digit, multiply by 10 and add next digit
    while(*psz != 0)
    {
        *this = *this * 10;
        nDigit = *psz - '0';
        *this = *this + nDigit;
        ++psz;
    }
}

// convert from integer
MPI::MPI(INT32 n)
{
    Clear();
    m_ar[0] = n;
}

/*****************************************************************************/
// ASSIGNMENTS
/*****************************************************************************/

MPI MPI::operator=(const MPI& m)
{
    // protect against assigning object to itself
    if(&m == this) return *this;

    // duplicate the array
    int i;
    for(i=0; i<MAX_ARRAY; ++i)
        m_ar[i] = m.m_ar[i];

    m_bOverflow = m.m_bOverflow;
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

// normal add MPI with MPI
// Algorithm 14.7
MPI MPI::operator+(const MPI& m) const
{
    MPI w; // result
    int c;    // carry
    int i;    // index

    for(i=0, c=0; i<MAX_ARRAY; ++i)
    {
        // add current digits
        w.m_ar[i] = m_ar[i] + m.m_ar[i] + c;

        // test for digit overflow
        if(w.m_ar[i] & MOD_VALUE)
        {
            w.m_ar[i] -= MOD_VALUE;
            c = 1;
        }
        else
            c = 0;
    }

    // if still carry, whole result is overflow
    w.m_bOverflow |= (c != 0);

    return w;
}

// normal add MPI with int
MPI MPI::operator+(const int n) const
{
    MPI w;
    int c;     // carry
    int i;     // index

    // add argument digit
    w.m_ar[0] = m_ar[0] + n;

    // digit overflow?
    c = 0;
    if(w.m_ar[0] & MOD_VALUE)                         
    {
        w.m_ar[0] -= MOD_VALUE;
        c = 1;
    }

    // copy all other digits, propagating possible carry
    for(i=1; i<MAX_ARRAY; ++i)
    {
        // copy digit, add carry
        w.m_ar[i] = m_ar[i] + c;

        // test for digit overflow
        if(w.m_ar[i] & MOD_VALUE)
        {
            w.m_ar[i] -= MOD_VALUE;
            c = 1;
        }
        else
            c = 0;
    }

    // if still carry, whole result is overflow
    w.m_bOverflow |= (c != 0);

    return w;
}

// normal subtract MPI with MPI
// Algorithm 14.9
MPI MPI::operator-(const MPI& m) const
{
    MPI w;
    int c;  // carry
    int i;  // index

    for(i=0, c=0; i<MAX_ARRAY; ++i)
    {
        // subtract current digits
        w.m_ar[i] = m_ar[i] - m.m_ar[i] - c;

        // test for digit overflow
        if(w.m_ar[i] & MOD_VALUE)
        {
            w.m_ar[i] += MOD_VALUE;
            c = 1;
        }
        else
            c = 0;
    }

    // if still carry, whole result is overflow
    w.m_bOverflow |= (c != 0);

    return w;
}

// normal subtract MPI with int
MPI MPI::operator-(const int n) const
{
    MPI w;
    int c;  // carry
    int i;  // index

    // subtract argument digit
    w.m_ar[0] = m_ar[0] - n;

    // digit overflow?
    c = 0;
    if(w.m_ar[0] & MOD_VALUE)                         
    {
        w.m_ar[0] += MOD_VALUE;
        c = 1;
    }

    for(i=1; i<MAX_ARRAY; ++i)
    {
        // subtract current digits
        w.m_ar[i] = m_ar[i] - c;

        // test for digit overflow
        if(w.m_ar[i] & MOD_VALUE)
        {
            w.m_ar[i] += MOD_VALUE;
            c = 1;
        }
        else
            c = 0;
    }

    // if still carry, whole result is overflow
    w.m_bOverflow |= (c != 0);

    return w;
}

/*****************************************************************************/
// MULTIPLICATION
/*****************************************************************************/

// normal multiply MPI with MPI
// Algorithm 14.12
//MPI MPI::MultALR(MPI& y) const
MPI MPI::operator*(const MPI& y) const
{
    MPI w;  // result
    int i;  // index
    int j;  // index
    int n;  // # digits in x
    int t;  // # digits in y

    INT64 uv; // double the precision of single digit
    INT64 c;  // carry

    n = Size();   // get # digits in x
    t = y.Size(); // get # digits in y

    for(i=0; i<t; ++i)
    {
        c = 0;

        for(j=0; j<n; ++j)
        {                    
            if(i+j > MAX_ARRAY-1)
                continue;

            uv = w.m_ar[i+j] + m_ar[j] * y.m_ar[i] + c;
            w.m_ar[i+j] = (INT32)(uv & (MOD_VALUE-1));  // LOWORD
            c = ((uv >> SHIFT_VALUE) & (MOD_VALUE-1));  // HIWORD 
        }

        if(i+j < MAX_ARRAY)
            w.m_ar[i+j] += (INT32)c;

        else if(c > 0)     
            w.m_bOverflow = TRUE;
    }
    return w;
}

// normal multiply MPI with int
MPI MPI::operator*(const int y) const
{
    MPI w;  // result
    int j;  // index
    int n;  // # digits in x

    INT64 uv; // double the precision of single digit
    INT64 c;  // carry

    n = Size()+1;   // get # digits in x + 1
    c = 0;

    for(j=0; j<n; ++j)
    {
        uv = (INT64)m_ar[j] * (INT64)y + c;
        w.m_ar[j] = (INT32)(uv & (MOD_VALUE-1));   // LOWORD
        c = ((uv >> SHIFT_VALUE) & (MOD_VALUE-1)); // HIWORD
    }

    if(c > 0)
        w.m_bOverflow = TRUE;

    return w;
}

// normal multiply MPI with MPI
// a la russe algorithm, brassard p. 4
//MPI MPI::operator*(const MPI& m) const
MPI MPI::MultALR(MPI& m) const
{
    MPI mProd; // result
    MPI x; // multiplicand
    MPI y; // multiplier
    int i, j;  // indices

    // optimize order, better if multiplier is smaller
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

    // calc how many loops to make
    j = x.Size() * SHIFT_VALUE;

    for(i=0; i<j; ++i)
    {
        // sum if odd
        if(y.m_ar[0] & 1)
            mProd += x;

        x.Mult2();
        y.Div2();
    }
    return mProd;
}

// normal multiply MPI with MPI
// divide and conquer algorithm, brassard p. 219-223
MPI MPI::MultDC(MPI& m) const
{
    MPI mProd;       // result
    MPI w, x, y, z;  // split pieces of the arguments
    MPI p, q, r;     // intermediate products
    int i;           // index
    int nFull, nHalf; // full size and split position of arguements

    nFull = Largest(m)-1;  // use helper fn()
    nHalf = nFull/2;

    // assign w
    for(i=nHalf; i<nFull; ++i)
        w.m_ar[i-nHalf] = m_ar[i];

    // assign x
    for(i=0; i<nHalf; ++i)
        x.m_ar[i] = m_ar[i];

    // assign y
    for(i=nHalf; i<nFull; ++i)
        y.m_ar[i-nHalf] = m.m_ar[i];

    // assign z
    for(i=0; i<nHalf; ++i)
        z.m_ar[i] = m.m_ar[i];

    // calculate three products
    p = w * y;
    q = x * z;
    r = (w + x) * (y + z);

    mProd = p;
    mProd.ShiftLeft(nHalf);

    mProd += r - p - q;
    mProd.ShiftLeft(nHalf);

    mProd += q;

    return mProd;
}

/*****************************************************************************/
// DIVISION AND MODULUS
/*****************************************************************************/

// normal divide MPI with MPI
// Algorithm D by Donald Knuth, p. 257
MPI MPI::operator/(const MPI& m) const
{
    MPI u;    // dividend
    MPI v;    // divisor
    MPI q;    // quotient
	MPI s;    // trial subtract amount
    INT32 qh; // trial quotient

    int j;    // index
    int n;    // # digits in u
    int t;    // # digits in v

    // work with local copies
    u = *this;
    v = m;

    // normalize
    while(v.MSDigit() < MOD_VALUE/2)
    {
        u.Mult2();
        v.Mult2();
    }

    // get # digits
    t = v.Size();
    n = u.Size();

    // main calculation loop
    for(j=n; j>=t; --j)
    {
        // calculate trial quotient, get LOWORD
		INT64 temp = (((INT64)u.m_ar[j] * MOD_VALUE) + u.m_ar[j-1]) / v.m_ar[t-1];
		qh = (INT32)temp;

        // adjust quotient if too large
        if(qh & MOD_VALUE)
            qh = MOD_VALUE -1;

        // multiply and subtract
        s = v;
        s.ShiftLeft(j-t);
		u -= s * qh;

		// adjust result if it went negative
        while(u.m_ar[MAX_ARRAY-1] == MOD_VALUE-1)
        {
            --qh;   // adjust quotient digit
            u += s; // add back
            u.m_bOverflow = FALSE;
        }

        // set the quotient digit we just found
        q.m_ar[j-t] = qh;
    }
    return q;
}

// normal divide MPI with int
MPI MPI::operator/(const int n) const
{
    return *this / MPI(n);
}


// normal modulus, MPI with MPI
// Algorithm D by Donald Knuth, p. 257
MPI MPI::operator%(const MPI& m) const
{
    MPI u;    // dividend
    MPI v;    // divisor
	MPI s;    // trial subtract amount
    INT32 qh; // trial quotient
    int d;    // normalization factor

    int j;    // index
    int n;    // # digits in u
    int t;    // # digits in v

    // work with local copies
    u = *this;
    v = m;

    // normalize
    d = 0;
    while(v.MSDigit() < MOD_VALUE/2)
    {
        u.Mult2();
        v.Mult2();
        ++d;
    }

    // get # digits
    t = v.Size();
    n = u.Size();

    // main calculation loop
    for(j=n; j>=t; --j)
    {
        // calculate trial quotient, get LOWORD
		INT64 temp = (((INT64)u.m_ar[j] * MOD_VALUE) + u.m_ar[j-1]) / v.m_ar[t-1];
		qh = (INT32)temp;

        // adjust quotient if too large
        if(qh & MOD_VALUE)
            qh = MOD_VALUE -1;

        // multiply and subtract
		s = v;
		s.ShiftLeft(j-t);
		u -= s * qh;

		// adjust result if it went negative
        while(u.m_ar[MAX_ARRAY-1] == MOD_VALUE-1)
            u += s; // add back
    }

    // un-normalize
    for(j=0; j<d; ++j)
        u.Div2();

	// fix flag if it ever went negative
	u.m_bOverflow = FALSE;

    return u;  // the remainder
}

// normal modulus MPI with int
MPI MPI::operator%(const int n) const
{
    return *this % MPI(n);
}

/*****************************************************************************/
// EXPONENTS AND MODULAR ARITHMETIC
/*****************************************************************************/

// exponential
MPI MPI::operator^(const MPI& y) const
{
    MPI w;  // return value
    MPI s;  // shifted exponent
    int n;  // digits in exponent
    int k;  // bits in exponent
    int i;  // index

    // shifted exponent
    s = y;

    // how many times to loop
    n = y.Size();
    k = n * SHIFT_VALUE;

    // x ^ 0 = 1
    if(!n) return MPI(1);

    // main loop
    for(i=0, w=1; i<k; ++i)
    {
        // square
        w *= w;

        // multiply
        if(s.m_ar[n-1] & MOD_VALUE/2)
            w *= *this;

        // shift
        s.Mult2();
    }
    return w;
}

// exponential
MPI MPI::operator^(const int n) const
{
    return *this ^ MPI(n);
}

// modular multiplication
// Algorihthm using definition of modular arithmetic
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

// modular exponetial
// Algorihthm using repeated squaring, CLR p. 829
MPI MPI::ModPow(const MPI& y, const MPI& m) const
{
    MPI w;  // return value
    MPI s;  // shifted exponent
    int n;  // digits in exponent
    int k;  // bits in exponent
    int i;  // index

    // shifted exponent
    s = y;

    // how many times to loop
    n = y.Size();
    k = n * SHIFT_VALUE;

    // x ^ 0 = 1
    if(!n) return MPI(1);

    // main loop
    for(i=0, w=1; i<k; ++i)
    {
        // square
        w *= w;
        w %= m;

        // multiply
        if(s.m_ar[n-1] & MOD_VALUE/2)
        {
            w *= *this;
            w %= m;
        }

        // shift
        s.Mult2();
    }
    return w;
}

/*****************************************************************************/
// ARITHMETIC SHORTCUT FORMS
/*****************************************************************************/

MPI MPI::operator+=(const MPI& m) { return *this = *this + m; }
MPI MPI::operator-=(const MPI& m) { return *this = *this - m; }
MPI MPI::operator*=(const MPI& m) { return *this = *this * m; }
MPI MPI::operator/=(const MPI& m) { return *this = *this / m; }
MPI MPI::operator%=(const MPI& m) { return *this = *this % m; }
MPI MPI::operator^=(const MPI& m) { return *this = *this ^ m; }

MPI MPI::operator+=(const int n)  { return *this = *this + n; }
MPI MPI::operator-=(const int n)  { return *this = *this - n; }
MPI MPI::operator*=(const int n)  { return *this = *this * n; }
MPI MPI::operator/=(const int n)  { return *this = *this / n; }
MPI MPI::operator%=(const int n)  { return *this = *this % n; }
MPI MPI::operator^=(const int n)  { return *this = *this ^ n; }

MPI MPI::operator++() { *this = *this + 1; return *this; }
MPI MPI::operator++(int) { MPI m = *this; *this = *this + 1; return m; }
MPI MPI::operator--() { *this = *this - 1; return *this; }
MPI MPI::operator--(int) { MPI m = *this; *this = *this - 1; return m; }

/*****************************************************************************/
// LOGICAL OPERATORS/ COMPARISON/ SHIFTING
/*****************************************************************************/

BOOL MPI::operator<(const MPI& m) const
{
    MPI x;

    x = *this - m;
    return x.m_bOverflow;
}

BOOL MPI::operator<=(const MPI& m) const
{
    MPI x;

    x = m - *this;
    return !x.m_bOverflow;
}

BOOL MPI::operator>(const MPI& m) const
{
    MPI x;

    x = m - *this;
    return x.m_bOverflow;
}

BOOL MPI::operator>=(const MPI& m) const
{
    MPI x;

    x = *this - m;
    return !x.m_bOverflow;
}

BOOL MPI::operator==(const MPI& m) const
{
    int i;

    for(i=0; i<MAX_ARRAY; ++i)
    {
        if(m.m_ar[i] != m_ar[i])
            return FALSE;
    }
    return TRUE;
}

BOOL MPI::operator!=(const MPI& m) const
{
    int i;

    for(i=0; i<MAX_ARRAY; ++i)
    {
        if(m.m_ar[i] != m_ar[i])
            return TRUE;
    }
    return FALSE;
}

/*****************************************************************************/
// DECIMAL CONVERSION AND I/O
/*****************************************************************************/

// convert the MPI type to a character string representation of decimal
char* MPI::String(char a[]) const
{
    int i;       // index into a[]
    MPI q;       // quotient
    MPI b;       // base
    MPI x;       // residue
    MPI zero;

    i = 0;
    b = 10;
    x = *this;

    // don't print junk
    if(Size() == MAX_ARRAY || m_bOverflow)
    {
        strcpy(a, "ERROR");
        return a;
    }

    do
    {
        q = x / b;
        x = (x - (q * b));
        a[i] = x.m_ar[0] + '0';
        ++i;
        x = q;
    }
    while(q != zero);

    a[i] = '\0';
    --i;

    // reverse the string into human readable order
    int j;
    char temp;
    for(j=0; j<=(i/2); ++j)
    {
        temp = a[j];
        a[j] = a[i-j];
        a[i-j] = temp;
    }
    return a;
}

// stream output
ostream& operator<<(ostream& os, MPI& m)
{
    char sz[MPI_BUFF];

    m.String(sz);
    os << sz;
    return os;
}

// stream input
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

// check for digit overflow and overflow flag
BOOL MPI::IsValid() const
{
    int i;

    if(m_bOverflow)
        return FALSE;

    for(i=0; i<MAX_ARRAY; ++i)
    {
        if(m_ar[i] > MOD_VALUE || m_ar[i] < 0)
            return FALSE;
    }

    return TRUE;
}

// diagnostic display internal representation
void MPI::Display() const
{
    int i;

    cout << "MPI [";
    for(i=0; i<MAX_ARRAY; ++i)
        cout << m_ar[i] << ", ";

    cout << "]\n";
}

// most significant digit
INT32 MPI::MSDigit() const
{
    int i;

    for(i=MAX_ARRAY-1; i>=0; --i)
    {
        if(m_ar[i] != 0)
            return m_ar[i];
    }
    return 0;
}

// how many significant digits
int MPI::Size() const
{
    int i;

    for(i=MAX_ARRAY-1; i>=0; --i)
    {
        if(m_ar[i] != 0)
            return i+1;
    }
    return 0;
}

// return the largest size of two arguments
inline int MPI::Largest(MPI& m) const
{
    int a, b, s;

    a = Size();   // size of left argument
    b = m.Size(); // size of right argument

    // find largest size argument
    if(a > b)
        s = a;
    else
        s = b;

    return s + 1;
}

// shift by n digit positions
void MPI::ShiftLeft(const int n)
{
    int i;  // index

    if(n >= MAX_ARRAY-1)
    {
        m_bOverflow = TRUE;
        return;
    }

    for(i=MAX_ARRAY-1; i>n-1; --i)
        m_ar[i] = m_ar[i-n];

    for(i=n-1; i>=0; --i)
        m_ar[i] = 0;
}

// shift by n digit positions,
void MPI::ShiftRight(const int n)
{   
    int i;  // index

    if(n >= MAX_ARRAY-1)
    {
        m_bOverflow = TRUE;
        return;
    }

    for(i=n; i<MAX_ARRAY; ++i)
        m_ar[i-n] = m_ar[i];

    for(i=MAX_ARRAY-n; i<MAX_ARRAY; ++i)
        m_ar[i] = 0;
}

// multiply using bit shift
void MPI::Mult2()
{
    int i;   // index
    INT32 c; // carry

    for(i=0, c=0; i<MAX_ARRAY; ++i)
    {
        m_ar[i] <<= 1;
        if(c) m_ar[i] |= 1;
        c = m_ar[i] & MOD_VALUE;
        m_ar[i] &= MOD_VALUE-1;
    }
}

// divide using bit shift
void MPI::Div2()
{
    int i;   // index
    INT32 c; // carry

    for(i=MAX_ARRAY-1, c=0; i>=0; --i)
    {
        if(c) m_ar[i] |= MOD_VALUE;
        c = m_ar[i] & 1;
        m_ar[i] >>= 1;
    }
}

