#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// return maximum value
#define MAX(i, j) (((i) > (j)) ? (i) : (j))

// 0-based index for square matrix
#define ZERO_INDEX(i, j, dim) ((i) * (dim) + (j))

// mymv : ret := A * x
// A : row-major order
void mymv(double ret[], double A[], double x[], int dim)
{
	int i, j;

	for(i = 0; i < dim; i++)
	{
		ret[i] = 0.0;
		for(j = 0; j < dim; j++)
			ret[i] += A[ZERO_INDEX(i, j, dim)] * x[j];
	}
}

// myaxpy : ret := alpha * x + y
void myaxpy(double ret[], double alpha, double x[], double y[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		ret[i] = alpha * x[i] + y[i];
}

// mydotp : x^T * y
double mydotp(double x[], double y[], int dim)
{
	int i;
	double ret = 0.0;

	for(i = 0; i < dim; i++)
		ret += x[i] * y[i];

	return ret;
}

// mynorm2 : ||x||_2
double mynorm2(double x[], int dim)
{
	int i;
	double ret = 0.0;

	for(i = 0; i < dim; i++)
		ret += x[i] * x[i];

	ret = sqrt(ret);

	return ret;
}

// mycopy : x := y
void mycopy(double x[], double y[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		x[i] = y[i];
}

// set0 : x := 0
void set0(double x[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		x[i] = 0.0;
}

// relerr: return relative error of val
double get_relerr(double val, double true_val)
{
	if(true_val != 0.0)
		return fabs((val - true_val) / true_val);
	else
		return fabs(val - true_val);
}

// Solve A * x = b
int conjugate_gradient(double x[], double A[], double b[], int dim, double rel_tol, double abs_tol, int maxitimes)
{
	int itimes;

	// vectors
	double *r, *x_k, *p, *w;

	// constant
	double alpha, beta, tmp_val, init_r_norm2, r_new_norm2;

	// Initialize
	r     = (double *)calloc(dim, sizeof(double));
	p     = (double *)calloc(dim, sizeof(double));
	w     = (double *)calloc(dim, sizeof(double));

	// x_0 := 0
	set0(x, dim);

	// r := b - A * x_0;
	mycopy(r, b, dim); // x_0 = 0

	// init_r_norm2 = ||r||_2
	init_r_norm2 = mynorm2(r, dim);

	// p := r
	mycopy(p, r, dim);

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// w := A * p
		mymv(w, A, p, dim);

		// alpha := (r, p) / (p, A * p)
		alpha = mydotp(r, p, dim);
		tmp_val = mydotp(p, w, dim);

		if(tmp_val == 0)
		{
			itimes = -1;
			break;
		}
		alpha /= tmp_val;

		// x := x + alpha * p
		myaxpy(x, alpha, p, x, dim);

		// r_new := r - alpha * A * p
		beta = mydotp(r, r, dim);
		myaxpy(r, -alpha, w, r, dim);

		// beta := ||r_new||_2^2 / ||r||_2^2
		r_new_norm2 = mynorm2(r, dim);
		beta = r_new_norm2 * r_new_norm2 / beta;

		// check residual
		if(r_new_norm2 <= rel_tol * init_r_norm2 + abs_tol)
			break;
		
		printf("%3d %10.3e %10.3e\n", itimes, alpha, r_new_norm2 / init_r_norm2);
	
		// p := r + beta * p
		myaxpy(p, beta, p, r, dim);

	}

	// clean
	free(r);
	free(p);
	free(w);

	return itimes;
}

// set matrix, true_x and b
void set_test_linear_eq(double A[], double true_x[], double b[], int dim)
{
	int i, j;

	// set A
	for(i = 0; i < dim; i++)
	{
		true_x[i] = (double)(i + 1);
		for(j = 0; j < dim; j++)
			A[ZERO_INDEX(i, j, dim)] = (double)dim - (double)MAX(i, j);
	}

	// b := A * true_x
	mymv(b, A, true_x, dim);

	// print the linear equation
	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
			printf("%10.3e ", A[ZERO_INDEX(i, j, dim)]);
		
		printf("  x  = %10.3e\n", b[i]);
	}
}

int main(int argc, char *argv[])
{
	int i, j, dimension, cg_itimes;
	double *matrix, *true_x, *b, *x;

	if(argc <= 1)
	{
		fprintf(stderr, "USAGE: %s [dimension]\n", argv[0]);
		return EXIT_SUCCESS;
	}

	dimension = atoi(argv[1]);

	if(dimension <= 1)
	{
		fprintf(stderr, "ERROR: dimension = %d is illegal!\n", dimension);
		return EXIT_FAILURE;
	}

	// initialize
	matrix = (double *)calloc(dimension * dimension, sizeof(double));
	true_x = (double *)calloc(dimension, sizeof(double));
	x      = (double *)calloc(dimension, sizeof(double));
	b      = (double *)calloc(dimension, sizeof(double));

	// set test problem
	set_test_linear_eq(matrix, true_x, b, dimension);

	// run conjugate-gradient routine
	cg_itimes = conjugate_gradient(x, matrix, b, dimension, 1.0e-10, 1.0e-100, dimension * 5);

	// print solution
	printf("cg iterative times: %d\n", cg_itimes);
	for(i = 0; i < dimension; i++)
		printf("%3d %25.17e %10.3e\n", i, x[i], get_relerr(x[i], true_x[i]));

	// free
	free(matrix);
	free(true_x);
	free(x);
	free(b);

	return EXIT_SUCCESS;
}
