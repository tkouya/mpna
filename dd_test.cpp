//******************************************************************************
// dd_test.cpp : Double-double precision computation with QD library
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

#define QD_INLINE
#include "qd/qd_real.h"
#include "qd/fpu.h"

using namespace std;

int main()
{
	int i, j;
	dd_real dd_a, dd_b;
	unsigned int old_cw;

	// DD
	fpu_fix_start(&old_cw);

	dd_a = sqrt((dd_real)2) * 2;
	dd_b = sqrt((dd_real)2);

	printf("dd_a[0] = %25.17e\n", dd_a.x[0]);
	printf("dd_a[1] = %25.17e\n", dd_a.x[1]);

	cout << "DD decimal digigs = " << dd_real::_ndigits << endl;

	cout << setprecision(dd_real::_ndigits) << "a     = " << dd_a << endl;
	cout << setprecision(dd_real::_ndigits) << "b     = " << dd_b << endl;

	cout << setprecision(dd_real::_ndigits) << "a + b = " << dd_a + dd_b << endl;
	cout << setprecision(dd_real::_ndigits) << "a - b = " << dd_a - dd_b << endl;
	cout << setprecision(dd_real::_ndigits) << "a * b = " << dd_a * dd_b << endl;
	cout << setprecision(dd_real::_ndigits) << "a / b = " << dd_a / dd_b << endl;

	fpu_fix_end(&old_cw);

	return 0;

}
