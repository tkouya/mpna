/**********************************************************/
/* Power Method with mpreal/MPFR/GMP                      */
/*                                                        */
/* Copyright (c) 2016 Tomonori Kouya, All rights reserved */
/* Version 0.0: 2016-12-20 (Tue) First published          */
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
//#include "template_linear.h"
#include "template_linear_power.h"

// Time routines
#include "get_secv.h"

int main(int argc, char *argv[])
{
	unsigned long prec;
	int i, j, dimension, power_itimes;
	mpreal *matrix, *x, min_eig, max_eig;
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
	x      = new mpreal[dimension];

	// set test problem
	//set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);
	set_test_linear_eq_hilbert<mpreal>(matrix, NULL, NULL, dimension);

	// run conjugate-gradient routine
	//rel_tol = "1.0e-10";
	rel_tol = "1.0e-100"; // dim = 128
	abs_tol = "1.0e-200";

	start_time = get_secv();
	power_itimes = power_method<mpreal>(max_eig, x, matrix, dimension, rel_tol, abs_tol, dimension * 5);
	end_time = get_secv() - start_time;

	// print solution
	cout << "-- mpreal(" << prec << " bits) precision --" << endl;
	cout << "dimension of prob.: " << dimension << endl;
	cout << "power iterative times: " << power_itimes << endl;
	cout << "power comp. time(sec): " << end_time << endl;
	cout << "max_eig : " << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << max_eig << endl;
	//for(i = 0; i < dimension; i++)
	//	cout << setw(3) << i << " " << scientific << setprecision(17) << x[i] << endl;

	start_time = get_secv();
	power_itimes = inverse_power_method_lu<mpreal>(min_eig, x, matrix, dimension, rel_tol, abs_tol, dimension * 5);
	end_time = get_secv() - start_time;

	// print solution
	cout << "-- mpreal(" << prec << " bits) precision --" << endl;
	cout << "dimension of prob.: " << dimension << endl;
	cout << "inverse power iterative times: " << power_itimes << endl;
	cout << "inverse power comp. time(sec): " << end_time << endl;
	cout << "min_eig : " << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << min_eig << endl;
	//for(i = 0; i < dimension; i++)
	//	cout << setw(3) << i << " " << scientific << setprecision(17) << x[i] << endl;

	// print condition number
	cout << "max_eig / min_eig = " << setprecision(mpfr::bits2digits(mpreal::get_default_prec())) << max_eig / min_eig << endl;

	// free
	delete_array<mpreal>(matrix, dimension * dimension);
	delete_array<mpreal>(x, dimension);

	return EXIT_SUCCESS;
}
