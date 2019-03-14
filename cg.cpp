#include <iostream>
#include <iomanip>

#include <cstdlib>
#include <cmath>

using namespace std;

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
	r = new double[dim];
	p = new double[dim];
	w = new double[dim];

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

		cout << setw(3) << itimes << " " << scientific << setprecision(3) << alpha << " " << r_new_norm2 / init_r_norm2 << endl;

		// p := r + beta * p
		myaxpy(p, beta, p, r, dim);

	}

	// clean
	delete r;
	delete p;
	delete w;

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
			cout << scientific << setprecision(3) << A[ZERO_INDEX(i, j, dim)] << " ";
		
		cout << "  x   " << scientific << setprecision(3) << b[i] << endl;
	}
}

int main(int argc, char *argv[])
{
	int i, j, dimension, cg_itimes;
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

	// set test problem
	set_test_linear_eq(matrix, true_x, b, dimension);

	// run conjugate-gradient routine
	cg_itimes = conjugate_gradient(x, matrix, b, dimension, 1.0e-10, 1.0e-100, dimension * 5);

	// print solution
	cout << "cg iterative times: " << cg_itimes << endl;
	for(i = 0; i < dimension; i++)
		cout << setw(3) << i << " " << scientific << setprecision(17) << x[i] << " " << setprecision(3) << get_relerr(x[i], true_x[i]) << endl;

	// free
	delete matrix;
	delete true_x;
	delete x;
	delete b;

	return EXIT_SUCCESS;
}
