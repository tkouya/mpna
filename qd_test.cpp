/********************************************************************************/
/* qd_test.cpp: Quad-double precision Linear Computation Library                */
/* Copyright (C) 2015-2016 Tomonori Kouya                                       */
/*                                                                              */
/* This program is free software: you can redistribute it and/or modify it      */
/* under the terms of the GNU Lesser General Public License as published by the */
/* Free Software Foundation, either version 3 of the License or any later       */
/* version.                                                                     */
/*                                                                              */
/* This program is distributed in the hope that it will be useful, but WITHOUT  */
/* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        */
/* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License */
/* for more details.                                                            */
/*                                                                              */
/* You should have received a copy of the GNU Lesser General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.        */
/*                                                                              */
/********************************************************************************/
#include <iostream>
#include <iomanip>

#define QD_INLINE
#include "qd/qd_real.h"
#include "qd/fpu.h"

using namespace std;

int main()
{
	int i, j;
	qd_real qd_a, qd_b;
	unsigned int old_cw;

	// QD
	fpu_fix_start(&old_cw);

	qd_a = sqrt((qd_real)2) * 2;
	qd_b = sqrt((qd_real)2);

	cout << "QD decimal digits = " << qd_real::_ndigits << endl;

	cout << setprecision(qd_real::_ndigits) << "a     = " << qd_a << endl;
	cout << setprecision(qd_real::_ndigits) << "b     = " << qd_b << "\n";

	cout << setprecision(qd_real::_ndigits) << "a + b = " << qd_a + qd_b << "\n";
	cout << setprecision(qd_real::_ndigits) << "a - b = " << qd_a - qd_b << "\n";
	cout << setprecision(qd_real::_ndigits) << "a * b = " << qd_a * qd_b << "\n";
	cout << setprecision(qd_real::_ndigits) << "a / b = " << qd_a / qd_b << "\n";

	fpu_fix_end(&old_cw);

	return 0;
}
