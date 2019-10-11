//******************************************************************************
// test_iterative_ref_mpreal.cpp : Interative refinement method based 
//                         on direct method (double, dd_real and qd_real + MPFR)
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

#include <cstdlib>
#include <cmath>

// Multiple precision with QD
#define QD_INLINE
#include "qd/qd_real.h"
#include "qd/fpu.h"

// mpreal
#include "mpreal.h"

// Template linear compucation with double, QD, MPFR/GMP
#include "template_linear.h"

// Time routines
#include "get_secv.h"
/* Simple Estimation of Condition Number */

using namespace std;
using namespace mpfr;

int main(int argc, char *argv[])
{
	unsigned long prec;
	int i, j, dimension, *pivot, itimes;
	mpreal *matrix, *true_x, *b, *x;
	mpreal rtol, atol;
	double start_time, end_time;
	unsigned int old_cw;

	if(argc <= 2)
	{
		cerr << "USAGE: " << argv[0] << " [dimension] [prec]" << endl;
		return EXIT_SUCCESS;
	}

	prec = atoi(argv[2]);
	if(prec <= 1)
	{
		cerr << "ERROR: prec = " << prec << " is illegal!" << endl;
		return EXIT_FAILURE;
	}
	mpreal::set_default_prec(prec);

	dimension = atoi(argv[1]);

	if(dimension <= 1)
	{
		cerr << "ERROR: dimension = " << dimension << " is illegal!" << endl;
		return EXIT_FAILURE;
	}

	fpu_fix_start(&old_cw);

	// initialize
	matrix = new mpreal[dimension * dimension];
	true_x = new mpreal[dimension];
	x      = new mpreal[dimension];
	b      = new mpreal[dimension];
	pivot  = new int[dimension];

	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);

	// run LU decomposion
	start_time = get_secv();
	LU<mpreal>(matrix, dimension, pivot);

	// backward & forward substitution
	solve_LU_linear_eq<mpreal>(x, matrix, b, dimension, pivot);
	end_time = get_secv() - start_time;

	cout << "comp.time(second): " << end_time << ", relerr = " << get_relerr_norm2(x, true_x, dimension) << endl;

	// double-MPFR
	// set test problem
	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);
//	rtol = "1.0e-300"; atol = "0.0";
	rtol = "1.0e-1000"; atol = "0.0";
//	rtol = "1.0e-50"; atol = "0.0";

	start_time = get_secv();
	itimes = iterative_refinement<mpreal, double>(x, matrix, b, rtol, atol, dimension, dimension * 10);
	end_time = get_secv() - start_time;
	cout << "Iterative Times: " << itimes << ", Dimension: " << dimension << endl;
	cout << "comp.time(second): " << end_time << ", relerr = " << get_relerr_norm2(x, true_x, dimension) << endl;

	// print solution
//	for(i = 0; i < dimension; i++)
//		cout << setw(3) << i << " " << scientific << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// DD-MPFR
	// set test problem
	start_time = get_secv();
	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);
//	rtol = "1.0e-100"; atol = "0.0";

	itimes = iterative_refinement<mpreal, dd_real>(x, matrix, b, rtol, atol, dimension, dimension * 10);
	end_time = get_secv() - start_time;
	cout << "Iterative Times: " << itimes << ", Dimension: " << dimension << endl;
	cout << "comp.time(second): " << end_time << ", relerr = " << get_relerr_norm2(x, true_x, dimension) << endl;

	// print solution
//	for(i = 0; i < dimension; i++)
//		cout << setw(3) << i << " " << scientific << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// QD-MPFR
	// set test problem
	start_time = get_secv();
	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);
//	rtol = "1.0e-200"; atol = "0.0";

	itimes = iterative_refinement<mpreal, qd_real>(x, matrix, b, rtol, atol, dimension, dimension * 10);
	end_time = get_secv() - start_time;
	cout << "Iterative Times: " << itimes << ", Dimension: " << dimension << endl;
	cout << "comp.time(second): " << end_time << ", relerr = " << get_relerr_norm2(x, true_x, dimension) << endl;

	// print solution
//	for(i = 0; i < dimension; i++)
//		cout << setw(3) << i << " " << scientific << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	fpu_fix_end(&old_cw);

	// free
	delete_array<mpreal>(matrix, dimension * dimension);
	delete_array<mpreal>(true_x, dimension);
	delete_array<mpreal>(x, dimension);
	delete_array<mpreal>(b, dimension);
	delete pivot;

	return EXIT_SUCCESS;
}
