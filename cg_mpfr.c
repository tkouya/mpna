//******************************************************************************
// cg_mpfr.c : Conjugate-gradient method with MPFR/GMP
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
#include "gmp.h"
#include "mpfr.h"

// C linear compucation with double, QD, MPFR/GMP
#include "linear_c.h"

int main(int argc, char *argv[])
{
	unsigned long prec;
	int i, j, dimension, cg_itimes;
	mpfr_t *matrix, *true_x, *b, *x;
	mpfr_t rel_tol, abs_tol, tmp_val;

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

	// initialize
	mpfr_inits2(prec, rel_tol, abs_tol, tmp_val, (mpfr_ptr)NULL);

	matrix = (mpfr_t *)calloc(dimension * dimension, sizeof(mpfr_t));
	true_x = (mpfr_t *)calloc(dimension, sizeof(mpfr_t));
	x      = (mpfr_t *)calloc(dimension, sizeof(mpfr_t));
	b      = (mpfr_t *)calloc(dimension, sizeof(mpfr_t));

	printf("callocs, dim = %d, prec = %ld\n", dimension, prec);

	mpfr_init2_array(matrix, dimension * dimension, prec);
	mpfr_init2_array(true_x, dimension, prec);
	mpfr_init2_array(x     , dimension, prec);
	mpfr_init2_array(b     , dimension, prec);

	// set test problem
	set_test_mpfr_linear_eq(matrix, true_x, b, dimension);

	printf("set_test_mpfr_linear_eq, dim = %d, prec = %ld\n", dimension, prec);

	// run conjugate-gradient routine
	mpfr_set_str(rel_tol, "1.0e-10", 10, _tk_default_rmode);
	mpfr_set_str(abs_tol, "1.0e-100", 10, _tk_default_rmode);

	cg_itimes = mpfr_conjugate_gradient(x, matrix, b, dimension, rel_tol, abs_tol, dimension * 5);

	// print solution
	printf("cg iterative times: %d\n", cg_itimes);
/*	for(i = 0; i < dimension; i++)
	{
		get_mpfr_relerr(tmp_val, x[i], true_x[i]);
		mpfr_printf("%3d %25.17RNe %10.3RNe\n", i, x[i], tmp_val);
	}
*/
	// free
	mpfr_clears(rel_tol, abs_tol, tmp_val, (mpfr_ptr)NULL);

	mpfr_clear_array(matrix, dimension * dimension);
	mpfr_clear_array(true_x, dimension);
	mpfr_clear_array(x, dimension);
	mpfr_clear_array(b, dimension);

	free(matrix);
	free(true_x);
	free(x);
	free(b);

	return EXIT_SUCCESS;
}
