//******************************************************************************
// template_lu_dd.cpp : Direct method with QD
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

// Multiple precision with QD
#define QD_INLINE
#include "qd/dd_real.h"
#include "qd/fpu.h"

// Template linear compucation using double, mpreal
#include "template_linear.h"

int main(int argc, char *argv[])
{
	unsigned int old_cw;
	int i, j, dimension, *pivot;
	dd_real *matrix, *true_x, *b, *x;
	dd_real rel_tol, abs_tol;

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

	set0<dd_real>(x, dimension);
	set0<dd_real>(b, dimension);

	// set test problem
	set_test_linear_eq<dd_real>(matrix, true_x, b, dimension);
	// QD must do with fpu_fix_start!
	fpu_fix_start(&old_cw);

	// run LU decomposion
	LU<dd_real>(matrix, dimension, pivot);

	// backward & forward substitution
	solve_LU_linear_eq<dd_real>(x, matrix, b, dimension, pivot);

	// print solution
	for(i = 0; i < dimension; i++)
		cout << setw(3) << i << " " << scientific << setprecision(dd_real::_ndigits) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// free
	delete_array<dd_real>(matrix, dimension * dimension);
	delete_array<dd_real>(true_x, dimension);
	delete_array<dd_real>(x, dimension);
	delete_array<dd_real>(b, dimension);
	delete pivot;

	// QD must end with fpu_fix_end!
	//fpu_fix_end(&old_cw);

	return EXIT_SUCCESS;
}
