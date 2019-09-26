//******************************************************************************
// complex_dd.cpp : Test program of double-double complex
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

#include "qd/qd_real.h"

using namespace std;

int main(void)
{
	complex<dd_real> a, b, c;

	// a = sqrt(2) + sqrt(3) * I
	// b = -sqrt(5) + pi * I

	a = complex<dd_real>(sqrt((dd_real)2), sqrt((dd_real)3));
	b = complex<dd_real>(-sqrt((dd_real)5), dd_real::_pi);

	cout << "a = " << a << endl;
	cout << "b = " << b << endl;

	cout << "a + b = " << a + b << endl;
	cout << "a - b = " << a - b << endl;
	cout << "a * b = " << a * b << endl;
	cout << "a / b = " << a / b << endl;
	cout << "sqrt(a) = " << sqrt(a) << endl;

	return 0;
}
