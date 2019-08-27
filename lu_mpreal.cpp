//******************************************************************************
// lu_dd.cpp : Direct method with MPFR C++ + MPFR/GMP
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

using namespace std;

// Multiple precision with MPFR/GMP
#include "mpreal.h"

using namespace mpfr;

// Template linear compucation using double, mpreal
#include "template_linear.h"

int main(int argc, char *argv[])
{
	unsigned long prec;
	int i, j, dimension, *pivot;
	mpreal *matrix, *true_x, *b, *x;
	mpreal rel_tol, abs_tol;

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

	// set test problem
	set_test_linear_eq<mpreal>(matrix, true_x, b, dimension);

	// run LU decomposion
	LU<mpreal>(matrix, dimension, pivot);

	// backward & forward substitution
	solve_LU_linear_eq<mpreal>(x, matrix, b, dimension, pivot);

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
