/**********************************************************/
/* Iterative Refinement Methods with various precision    */
/*                                                        */
/* Copyright (c) 2018 Tomonori Kouya, All rights reserved */
/*                                                        */
/* Version 0.0: 2018-12-15 (Tue) First published          */
/**********************************************************/
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
#include "template_linear_complete.h"

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

	// initialize
	matrix = new mpreal[dimension * dimension];
	true_x = new mpreal[dimension];
	x      = new mpreal[dimension];
	b      = new mpreal[dimension];
	pivot  = new int[dimension];

	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);

	// run LU decomposion
	//LU<mpreal>(matrix, dimension, pivot);

	// backward & forward substitution
	//solve_LU_linear_eq<mpreal>(x, matrix, b, dimension, pivot);

	// double-MPFR
	// set test problem
	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);
	rtol = "1.0e-50"; atol = "0.0";

	itimes = gesirsv<mpreal, double>(x, matrix, b, rtol, atol, dimension, dimension * 10);
	cout << "Iterative Times: " << itimes << ", Dimension: " << dimension << endl;

	// print solution
	for(i = 0; i < dimension; i++)
		cout << setw(3) << i << " " << scientific << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// DD-MPFR
	// set test problem
	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);
	rtol = "1.0e-100"; atol = "0.0";

	itimes = gesirsv<mpreal, dd_real>(x, matrix, b, rtol, atol, dimension, dimension * 10);
	cout << "Iterative Times: " << itimes << ", Dimension: " << dimension << endl;

	// print solution
	for(i = 0; i < dimension; i++)
		cout << setw(3) << i << " " << scientific << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// QD-MPFR
	// set test problem
	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);
	rtol = "1.0e-200"; atol = "0.0";

	itimes = gesirsv<mpreal, qd_real>(x, matrix, b, rtol, atol, dimension, dimension * 10);
	cout << "Iterative Times: " << itimes << ", Dimension: " << dimension << endl;

	// print solution
	for(i = 0; i < dimension; i++)
		cout << setw(3) << i << " " << scientific << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// free
	delete_array<mpreal>(matrix, dimension * dimension);
	delete_array<mpreal>(true_x, dimension);
	delete_array<mpreal>(x, dimension);
	delete_array<mpreal>(b, dimension);
	delete pivot;

	return EXIT_SUCCESS;
}
