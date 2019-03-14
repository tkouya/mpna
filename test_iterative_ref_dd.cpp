/**********************************************************/
/* Iterative Refinement Methods with DD precision         */
/*                                                        */
/* Copyright (c) 2018 Tomonori Kouya, All rights reserved */
/*                                                        */
/* Version 0.0: 2018-12-18 (Tue) First published          */
/**********************************************************/
#include <iostream>
#include <iomanip>

#include <cstdlib>
#include <cmath>

// Multiple precision with QD
#define QD_INLINE
#include "qd/qd_real.h"
#include "qd/fpu.h"

// Template linear compucation with double, QD, MPFR/GMP
#include "template_linear_complete.h"

// Time routines
#include "get_secv.h"
/* Simple Estimation of Condition Number */

using namespace std;

int main(int argc, char *argv[])
{
	unsigned long prec;
	int i, j, dimension, *pivot, itimes;
	dd_real *matrix, *true_x, *b, *x;
	dd_real rtol, atol;

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

	// initialize
	matrix = new dd_real[dimension * dimension];
	true_x = new dd_real[dimension];
	x      = new dd_real[dimension];
	b      = new dd_real[dimension];
	pivot  = new int[dimension];

	//set_test_linear_eq<qd_real>(matrix, true_x, b, dimension);

	// run LU decomposion
	//LU<mpreal>(matrix, dimension, pivot);

	// backward & forward substitution
	//solve_LU_linear_eq<mpreal>(x, matrix, b, dimension, pivot);

	// double-DD
	// set test problem
	set_test_linear_eq<dd_real>(matrix, true_x, b, dimension);
	rtol = "1.0e-50"; atol = "0.0";

	itimes = iterative_refinement<dd_real, double>(x, matrix, b, rtol, atol, dimension, dimension * 10);
	cout << "Iterative Times: " << itimes << ", Dimension: " << dimension << endl;

	// print solution
	for(i = 0; i < dimension; i++)
		cout << setw(3) << i << " " << scientific << setprecision(32) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// free
	delete_array<dd_real>(matrix, dimension * dimension);
	delete_array<dd_real>(true_x, dimension);
	delete_array<dd_real>(x, dimension);
	delete_array<dd_real>(b, dimension);
	delete pivot;

	return EXIT_SUCCESS;
}
