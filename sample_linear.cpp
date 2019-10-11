//******************************************************************************
// sample_linear.cpp : Sample code of multiple precision linear computation
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

// Template linear compucation with double, QD, MPFR/GMP
#include "template_linear.h"

int main(int argc, char *argv[])
{
	unsigned long prec;
	int i, j, dimension;
	mpreal *matrix, *b, *x;
	gmp_randstate_t state; // random

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

	// random number
	gmp_randinit_default(state);
	gmp_randseed_ui(state, (unsigned long)dimension); // seed := dimension

	// initialize
	matrix = new mpreal[dimension * dimension];
	x      = new mpreal[dimension];
	b      = new mpreal[dimension];

	// set x and matrix
	for(i = 0; i < dimension; i++)
	{
		x[i] = urandomb(state);
		for(j = 0; j < dimension; j++)
			matrix[ZERO_INDEX(i, j, dimension)] = urandomb(state);

	}

	// b := matrix * x
	mymv<mpreal>(b, matrix, x, dimension);

	// print the linear equation
	for(i = 0; i < dimension; i++)
	{
		for(j = 0; j < dimension; j++)
			cout << scientific << setprecision(17) << matrix[ZERO_INDEX(i, j, dimension)] << " ";
		
		cout << "  x   " << scientific << setprecision(17) << x[i] <<  " = " << b[i] << endl;
	}

	// free
	gmp_randclear(state);
	delete_array<mpreal>(matrix, dimension * dimension);
	delete_array<mpreal>(x, dimension);
	delete_array<mpreal>(b, dimension);

	return EXIT_SUCCESS;
}
