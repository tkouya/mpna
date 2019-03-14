/**********************************************************/
/* Conjugate-Gradient Linear Solver with mpreal/MPFR/GMP  */
/*                                                        */
/* Copyright (c) 2016 Tomonori Kouya, All rights reserved */
/* Version 0.0: 2016-11-17 (Thu) First published          */
/**********************************************************/
#include <iostream>
#include <iomanip>

#include <cstdlib>
#include <cmath>

using namespace std;

// Multiple precision with MPFR/GMP
#include "mpreal.h"

using namespace mpfr;

// Template linear compucation with double, QD, MPFR/GMP
#include "template_linear.h"

// Time routines
#include "get_secv.h"

int main(int argc, char *argv[])
{
	unsigned long prec;
	int i, j, dimension, cg_itimes;
	mpreal *matrix, *true_x, *b, *x;
	mpreal rel_tol, abs_tol;
	double start_time, end_time;

	if(argc <= 2)
	{
		cerr << "USAGE: " << argv[0] << " [dimension]  [prec]" << endl;
		return EXIT_SUCCESS;
	}

	// get precitsion in bits
	prec = (unsigned long)atoi(argv[2]);

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

	// set test problem
	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);

	// run conjugate-gradient routine
	rel_tol = "1.0e-10";
	abs_tol = "1.0e-200";

	start_time = get_secv();
	cg_itimes = conjugate_gradient<mpreal>(x, matrix, b, dimension, rel_tol, abs_tol, dimension * 5);
	end_time = get_secv() - start_time;

	// print solution
	cout << "-- mpreal(" << prec << " bits) precision --" << endl;
	cout << "dimension of prob.: " << dimension << endl;
	cout << "cg iterative times: " << cg_itimes << endl;
	cout << "cg comp. time(sec): " << end_time << endl;
	cout << "relerr_norm2      : " << get_relerr_norm2<mpreal>(x, true_x, dimension) << endl;
	//for(i = 0; i < dimension; i++)
	//	cout << setw(3) << i << " " << scientific << setprecision(17) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// free
	delete_array<mpreal>(matrix, dimension * dimension);
	delete_array<mpreal>(true_x, dimension);
	delete_array<mpreal>(x, dimension);
	delete_array<mpreal>(b, dimension);

	return EXIT_SUCCESS;
}
