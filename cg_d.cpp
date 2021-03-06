//******************************************************************************
// cg_d.cpp :  Conjugate-Gradient Linear Solver with double precision
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

// Template linear compucation with double, QD, MPFR/GMP
#include "template_linear.h"

// Time routines
#include "get_secv.h"

int main(int argc, char *argv[])
{
	int i, j, dimension, cg_itimes;
	double *matrix, *true_x, *b, *x;
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

	// initialize
	matrix = new double[dimension * dimension];
	true_x = new double[dimension];
	x      = new double[dimension];
	b      = new double[dimension];

	// set test problem
	set_test_linear_eq<double>(matrix, true_x, b, dimension);

	// run conjugate-gradient routine
	start_time = get_secv();
	cg_itimes = conjugate_gradient<double>(x, matrix, b, dimension, 1.0e-10, 1.0e-100, dimension * 5);
	end_time = get_secv() - start_time;

	// print solution
	cout << "-- double precision --" << endl;
	cout << "dimension of prob.: " << dimension << endl;
	cout << "cg iterative times: " << cg_itimes << endl;
	cout << "cg comp. time(sec): " << end_time << endl;
	cout << "relerr_norm2      : " << get_relerr_norm2<double>(x, true_x, dimension) << endl;
/*	for(i = 0; i < dimension; i++)
		cout << setw(3) << i << " " << scientific << setprecision(17) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;
*/
	// free
	delete_array<double>(matrix, dimension * dimension);
	delete_array<double>(true_x, dimension);
	delete_array<double>(x, dimension);
	delete_array<double>(b, dimension);

	return EXIT_SUCCESS;
}
