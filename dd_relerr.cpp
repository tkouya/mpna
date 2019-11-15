//******************************************************************************
// dd_relerr.cpp : Sample code to get relative error of dd_real variables
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

// Relative error of dd_real value
void dd_relerr(dd_real &relerr, dd_real approx, qd_real true_val)
{
	qd_real tmp_qd;

	tmp_qd = (qd_real)approx - true_val;
	if(true_val != 0)
		tmp_qd /= true_val;
	
	tmp_qd = abs(tmp_qd);
	relerr.x[0] = tmp_qd.x[0];
	relerr.x[1] = tmp_qd.x[1];
}

int main()
{
	int i, j;
	dd_real dd_a, dd_b, relerr;
	qd_real qd_a, qd_b;
	unsigned int old_cw;

	fpu_fix_start(&old_cw);

	// DD
	dd_a = sqrt((dd_real)2) * 2;
	dd_b = sqrt((dd_real)2);

	// QD
	qd_a = sqrt((qd_real)2) * 2;
	qd_b = sqrt((qd_real)2);

	cout <<  "                 dd_real value                Relative Error" << endl;

	dd_relerr(relerr, dd_a, qd_a);
	cout << "a     = " << setprecision(dd_real::_ndigits) << dd_a << " " << setprecision(3) << relerr << endl;
	dd_relerr(relerr, dd_b, qd_b);
	cout << "b     = " << setprecision(dd_real::_ndigits) << dd_b << " " << setprecision(3) << relerr << endl;

	dd_relerr(relerr, dd_a + dd_b, qd_a + qd_b);
	cout << setprecision(dd_real::_ndigits) << "a + b = " << dd_a + dd_b << " " << setprecision(3) << relerr << endl;

	dd_relerr(relerr, dd_a - dd_b, qd_a - qd_b);
	cout << setprecision(dd_real::_ndigits) << "a - b = " << dd_a - dd_b << " " << setprecision(3) << relerr << endl;

	dd_relerr(relerr, dd_a * dd_b, qd_a * qd_b);
	cout << setprecision(dd_real::_ndigits) << "a * b = " << dd_a * dd_b << " " << setprecision(3) << relerr << endl;

	dd_relerr(relerr, dd_a / dd_b, qd_a / qd_b);
	cout << setprecision(dd_real::_ndigits) << "a / b = " << dd_a / dd_b << " " << setprecision(3) << relerr << endl;

	fpu_fix_end(&old_cw);

}
