//******************************************************************************
// mpz_gcd_lcm.cpp : Sample code to get GCD and LCM of multiple precision 
// integers
// Copyright (C) 2019 Tomonori Kouya
// 
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License or any later
// version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
//******************************************************************************
#include <iostream>
#include "gmpxx.h"

using namespace std;

int main(void)
{
    mpz_class a, b, c;
    mpz_class s, t;

    a = "123456789012345678901234567890";
    b = 9876543210;

    // GCD(a, b)
    c = gcd(a, b);
    cout << "GCD(" << a << ", " << b << ") = " << c << endl;

    // LCM(a, b)
    c = lcm(a, b);
    cout << "LCM(" << a << ", " << b << ") = " << c << endl;

    // GCDEXT(a, b) -> a * s + b * t = c
    mpz_gcdext(c.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    cout << a << " * " << s << " + " << b << " * " << t << " = " << c << endl;

    return 0;
}
