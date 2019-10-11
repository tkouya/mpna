//******************************************************************************
// template_lu.cpp : Direct method with double precision
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

// Template linear compucation using double, mpreal
#include "template_linear.h"

int main(int argc, char *argv[])
{
	int i, j, dimension, cg_itimes, *pivot;
	double *matrix, *true_x, *b, *x;

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
	pivot  = new int[dimension];

	// set test problem
	set_test_linear_eq<double>(matrix, true_x, b, dimension);

	// run LU decomposion
	LU<double>(matrix, dimension, pivot);

	// backward & forward substitution
	solve_LU_linear_eq<double>(x, matrix, b, dimension, pivot);

	// print solution
	for(i = 0; i < dimension; i++)
		cout << setw(3) << i << " " << pivot[i] << " " << scientific << setprecision(17) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// free
	delete matrix;
	delete true_x;
	delete x;
	delete b;
	delete pivot;

	return EXIT_SUCCESS;
}
