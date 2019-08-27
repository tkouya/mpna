//******************************************************************************
// cg_double_omp.c : Conjugate-gradient method with double precision and OpenMP
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

// OpenMP enable ?
#ifdef _OPENMP
	// C linear compucation with double, QD, MPFR/GMP
	#include "linear_c_omp.h"
#else
	// C linear compucation with double, QD, MPFR/GMP
	#include "linear_c.h"
#endif // _OPENMP

// Time routines
#include "get_secv.h"

int main(int argc, char *argv[])
{
	int i, j, dimension, cg_itimes, num_threads;
	double *matrix, *true_x, *b, *x;
	double start_time, end_time;

	#ifdef _OPENMP
	printf("-- OpenMP enable! #procs = %d --\n", omp_get_num_procs());
	if(argc <= 2)
	{
		fprintf(stderr, "USAGE: %s [dimension] [#threads (<= %d) ]\n", argv[0], omp_get_num_procs());
		return EXIT_SUCCESS;
	}

	num_threads = atoi(argv[2]);
	omp_set_num_threads(num_threads);
	printf("#Threads = %d\n", omp_get_max_threads());

	#else // _OPENMP
	if(argc <= 1)
	{
		fprintf(stderr, "USAGE: %s [dimension]\n", argv[0]);
		return EXIT_SUCCESS;
	}
	#endif // _OPENMP

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

	// set test problem
	set_test_d_linear_eq(matrix, true_x, b, dimension);

	// run conjugate-gradient routine
	start_time = get_real_secv();
	cg_itimes = d_conjugate_gradient(x, matrix, b, dimension, 1.0e-10, 1.0e-100, dimension * 5);
	end_time = get_real_secv() - start_time;

	// print solution
	printf("-- double(53 bits) precision --\n");
	printf("dimension of prob.: %d\n", dimension);
	printf("cg iterative times: %d\n", cg_itimes);
	printf("cg comp. time(sec): %f\n", end_time);
	printf("relerr_norm2      : %10.3e\n", get_d_relerr_norm2(x, true_x, dimension));
	//for(i = 0; i < dimension; i++)
	//	printf("%3d %25.17e %10.3e\n", i, x[i], get_d_relerr(x[i], true_x[i]));

	// free
	free(matrix);
	free(true_x);
	free(x);
	free(b);

	return EXIT_SUCCESS;
}
