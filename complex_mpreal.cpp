//******************************************************************************
// complex_mpreal.cpp : Complex arithmetic with MPFR
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
#include <iomanip>
#include <complex>

#include "mpreal.h"

using namespace std;
using namespace mpfr;

int main(void)
{

	mpreal::set_default_prec(128);
	complex<mpreal> a, b, c;
	unsigned int dprec;

	// a = sqrt(2) + sqrt(3) * I
	// b = sqrt(3) - sqrt(5) * I

	dprec = (unsigned int)ceil(mpreal::get_default_prec() * log10(2.0));
	cout << "prec in bits: " << mpreal::get_default_prec() << endl;
	cout << "prec in dec.: " << dprec << endl;

	a = complex<mpreal>(sqrt((mpreal)2), sqrt((mpreal)3));
	b = complex<mpreal>(-sqrt((mpreal)5), const_pi());

	cout << "a = " << a << endl;
	cout << "b = " << b << endl;

	cout << "a + b = " << setprecision(dprec) << a + b << endl;
	cout << "a - b = " << setprecision(dprec) << a - b << endl;
	cout << "a * b = " << setprecision(dprec) << a * b << endl;
	cout << "a / b = " << setprecision(dprec) << a / b << endl;
	cout << "sqrt(a) = " << sqrt(a) << endl;

	return 0;
}