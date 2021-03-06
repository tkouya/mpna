//******************************************************************************
// lu_d.c : Direct method with double precision
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// C linear compucation with double, QD, MPFR/GMP
#include "linear_c.h"

int main(int argc, char *argv[])
{
	int i, j, dimension, *pivot;
	double *matrix, *true_x, *b, *x;

	if(argc <= 1)
	{
		fprintf(stderr, "USAGE: %s [dimension]\n", argv[0]);
		return EXIT_SUCCESS;
	}

	dimension = atoi(argv[1]);

	if(dimension <= 1)
	{
		fprintf(stderr, "ERROR: dimension = %d is illegal!", dimension);
		return EXIT_FAILURE;
	}

	// initialize
	matrix = (double *)calloc(dimension * dimension, sizeof(double));
	true_x = (double *)calloc(dimension, sizeof(double));
	x      = (double *)calloc(dimension, sizeof(double));
	b      = (double *)calloc(dimension, sizeof(double));
	pivot  = (int *)calloc(dimension, sizeof(int));

	// set test problem
	set_test_d_linear_eq(matrix, true_x, b, dimension);

	// LU decomposition & forward-backsubstitution
	d_LU(matrix, dimension, pivot);
	d_solve_LU_linear_eq(x, matrix, b, dimension, pivot);

	// print solution
	for(i = 0; i < dimension; i++)
		printf("%3d %25.17e %25.17e %10.3e\n", i, x[i], true_x[i], get_d_relerr(x[i], true_x[i]));

	// free
	free(matrix);
	free(true_x);
	free(x);
	free(b);
	free(pivot);

	return EXIT_SUCCESS;
}
