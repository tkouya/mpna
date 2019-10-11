//******************************************************************************
// sample_linear.c : Sample code of multiple precision linear computation
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

// MPFR/GMP is used
#include "mpfr.h"

// C linear compucation with double, MPFR/GMP
#include "linear_c.h"

int main(int argc, char *argv[])
{
	unsigned long prec;
	int i, j, dimension;
	mpfr_t *matrix, *b, *x;
	gmp_randstate_t state; // random

	if(argc <= 2)
	{
		fprintf(stderr, "USAGE: %s [dimension] [prec]\n", argv[0]);
		return EXIT_SUCCESS;
	}

	// get precitsion in bits
	prec = (unsigned long)atoi(argv[2]);
	dimension = atoi(argv[1]);

	if(dimension <= 1)
	{
		fprintf(stderr, "ERROR: dimension = %d is illegal!", dimension);
		return EXIT_FAILURE;
	}

	// random number
	gmp_randinit_default(state);
	gmp_randseed_ui(state, (unsigned long)dimension); // seed := dimension

	// initialize
	matrix = (mpfr_t *)calloc(dimension * dimension, sizeof(mpfr_t));
	x      = (mpfr_t *)calloc(dimension, sizeof(mpfr_t));
	b      = (mpfr_t *)calloc(dimension, sizeof(mpfr_t));

	printf("callocs, dim = %d, prec = %ld\n", dimension, prec);

	mpfr_init2_array(matrix, dimension * dimension, prec);
	mpfr_init2_array(x     , dimension, prec);
	mpfr_init2_array(b     , dimension, prec);

	// set x and matrix as random 
	for(i = 0; i < dimension; i++)
	{
		mpfr_urandomb(x[i], state);
		for(j = 0; j < dimension; j++)
			mpfr_urandomb(matrix[ZERO_INDEX(i, j, dimension)], state);
	}

	// b := matrix * _x
	mpfr_mymv(b, matrix, x, dimension);

	// print the linear equation
	for(i = 0; i < dimension; i++)
	{
		for(j = 0; j < dimension; j++)
			mpfr_printf("%15.7RNe ", matrix[ZERO_INDEX(i, j, dimension)]);
		
		mpfr_printf("  x   %15.7RNe = %15.7RNe\n", x[i], b[i]);
	}

	// free
	gmp_randclear(state);
	mpfr_clear_array(matrix, dimension * dimension);
	mpfr_clear_array(x, dimension);
	mpfr_clear_array(b, dimension);

	free(matrix);
	free(x);
	free(b);

	return EXIT_SUCCESS;
}
