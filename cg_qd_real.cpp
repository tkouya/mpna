/**********************************************************/
/* Conjugate-Gradient Linear Solver with qd_real          */
/*                                                        */
/* Copyright (c) 2016 Tomonori Kouya, All rights reserved */
/* Version 0.0: 2016-11-17 (Thu) First published          */
/**********************************************************/
#include <iostream>
#include <iomanip>

#include <cstdlib>
#include <cmath>

using namespace std;

// Multiple precision with QD
#include "qd/qd_real.h"

// Template linear compucation with double, QD, MPFR/GMP
#include "template_linear.h"

// Time routines
#include "get_secv.h"

int main(int argc, char *argv[])
{
	unsigned int old_cw;
	int i, j, dimension, cg_itimes;
	qd_real *matrix, *true_x, *b, *x;
	qd_real rel_tol, abs_tol;
	double start_time, end_time;

	if(argc <= 1)
	{
		cerr << "USAGE: " << argv[0] << " [dimension]" << endl;
		return EXIT_SUCCESS;
	}

	dimension = atoi(argv[1]);

	if(dimension <= 1)
	{
		cerr << "ERROR: dimension = " << dimension << " is illegal!" << endl;
		return EXIT_FAILURE;
	}

	// QD must do with fpu_fix_start!
	fpu_fix_start(&old_cw);

	// initialize
	matrix = new qd_real[dimension * dimension];
	true_x = new qd_real[dimension];
	x      = new qd_real[dimension];
	b      = new qd_real[dimension];

	// set test problem
	set_test_linear_eq<qd_real>(matrix, true_x, b, dimension);

	// run conjugate-gradient routine
	rel_tol = 1.0e-10;
	abs_tol = 1.0e-100;

	start_time = get_secv();
	cg_itimes = conjugate_gradient<qd_real>(x, matrix, b, dimension, rel_tol, abs_tol, dimension * 5);
	end_time = get_secv() - start_time;

	// print solution
	cout << "-- qd_real precision --" << endl;
	cout << "dimension of prob.: " << dimension << endl;
	cout << "cg iterative times: " << cg_itimes << endl;
	cout << "cg comp. time(sec): " << end_time << endl;
	cout << "relerr_norm2      : " << get_relerr_norm2<qd_real>(x, true_x, dimension) << endl;
	//for(i = 0; i < dimension; i++)
	//	cout << setw(3) << i << " " << scientific << setprecision(qd_real::_ndigits) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// free
	delete_array<qd_real>(matrix, dimension * dimension);
	delete_array<qd_real>(true_x, dimension);
	delete_array<qd_real>(x, dimension);
	delete_array<qd_real>(b, dimension);

	// QD must end with fpu_fix_end!
	fpu_fix_end(&old_cw);

	return EXIT_SUCCESS;
}
