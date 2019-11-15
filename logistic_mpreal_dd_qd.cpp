//******************************************************************************
// logistic_mpreal_dd_qd.cpp : Sequence calculation based on logistic map 
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

// QD
#include "qd_real.h"

// Multiple precision with MPFR/GMP
#include "mpreal.h"

using namespace std;
using namespace mpfr;

// mpfr_[set,get]_[dd,qd]
#include "mpfr_dd_qd.h"


int main(int argc, char *argv[])
{
	int i;
	int num_bits, num_decimal;

	// Check arguments
	if(argc <= 1)
	{
		cerr << "Usage: " << argv[0] << " [num_bits]" << endl;
		return 0;
	}

	// set default precision
	num_bits = atoi(argv[1]);
	if(num_bits <= 24)
		num_bits = 24;

	num_decimal = (int)ceil(log10(2.0) * (double)num_bits);
	mpreal::set_default_prec(num_bits);

	cout << "num_bits = " << num_bits << ", num_decimal = " << num_decimal << endl;

	dd_real dd_x[102], mpfr_ddreal; // dd
	qd_real qd_x[102], mpfr_qdreal; // qd
	mpreal x[102], dd_mpfr, qd_mpfr, reldiff_dd, reldiff_qd, reldiff_mpfr; // mpfr
	mpreal long_x[102]; // longer
	double mpfr_dd[2], mpfr_qd[4];

	// initialize as twice precision
	for(i = 0; i < 102; i++)
		long_x[i].set_prec(num_bits * 2);

	// fix FPU mode for QD
	fpu_fix_start(NULL);

	// set initial values
	x[0] = "0.7501";
	dd_x[0] = "0.7501";
	qd_x[0] = "0.7501";
	long_x[0] = "0.7501";

	for(i = 0; i <= 100; i++)
	{
		x[i + 1] = 4 * x[i] * (1 - x[i]);
		dd_x[i + 1] = 4 * dd_x[i] * (1 - dd_x[i]);
		qd_x[i + 1] = 4 * qd_x[i] * (1 - qd_x[i]);
		long_x[i + 1] = 4 * long_x[i] * (1 - long_x[i]);

		if((i % 10) == 0)
		{
			//cout << setw(5) << scientific << i << setprecision(num_decimal) << ", " << x[i] << setprecision(32) << ", " << dd_x[i] << endl;
			//cout << setw(5) << scientific << i << setprecision(num_decimal) << ", " << x[i] << setprecision(32) << ", " << dd_x[i] << ", " << qd_x[i] << endl;

			// dd := mpfr
			mpfr_get_dd(mpfr_dd, x[i].mpfr_srcptr(), MPFR_RNDN);
			mpfr_ddreal = dd_real(mpfr_dd[0], mpfr_dd[1]);
			//cout << setw(5) << scientific << i << setprecision(32) << ", " << mpfr_ddreal << setprecision(32) << ", " << dd_x[i] << endl;
			//cout << setw(5) << scientific << i << setprecision(32) << ", " << x[i] << endl;
			//cout << setw(5) << scientific << i << setprecision(32) << ", " << mpfr_ddreal << endl;

			// mpfr := dd
			mpfr_set_dd(dd_mpfr.mpfr_ptr(), dd_x[i].x, MPFR_RNDN);

			//cout << setw(5) << scientific << i << setprecision(32) << ", " << dd_mpfr << endl;
			//cout << setw(5) << scientific << i << setprecision(32) << ", " << dd_x[i] << endl;

			// qd := mpfr
			mpfr_get_qd(mpfr_qd, x[i].mpfr_srcptr(), MPFR_RNDN);
			mpfr_qdreal = qd_real(mpfr_qd[0], mpfr_qd[1], mpfr_qd[2], mpfr_qd[3]);

			//cout << setw(5) << scientific << i << setprecision(64) << ", " << x[i] << endl;
			//cout << setw(5) << scientific << i << setprecision(64) << ", " << mpfr_qdreal << endl;

			// mpfr := qd
			mpfr_set_qd(qd_mpfr.mpfr_ptr(), qd_x[i].x, MPFR_RNDN);

			//cout << setw(5) << scientific << i << setprecision(64) << ", " << qd_mpfr << endl;
			//cout << setw(5) << scientific << i << setprecision(64) << ", " << qd_x[i] << endl;

			// reldiff_xd := |x[i] - xd_mpfr| / x[i]
			mpfr_reldiff(reldiff_dd.mpfr_ptr(), long_x[i].mpfr_srcptr(), dd_mpfr.mpfr_srcptr(), MPFR_RNDN);
			mpfr_reldiff(reldiff_qd.mpfr_ptr(), long_x[i].mpfr_srcptr(), qd_mpfr.mpfr_srcptr(), MPFR_RNDN);
			mpfr_reldiff(reldiff_mpfr.mpfr_ptr(), long_x[i].mpfr_srcptr(), x[i].mpfr_srcptr(), MPFR_RNDN);

			cout << setw(5) << scientific << i << setprecision(5) << ", " << reldiff_mpfr << ", " << reldiff_dd << ", " << reldiff_qd << endl;
		}
	}

	return 0;
}
