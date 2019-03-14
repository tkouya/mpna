/**********************************************************/
/* Template for Multiple precision linear computation     */
/* Supporting double, dd_real, qd_real and mpreal         */
/*                                                        */
/* Copyright (c) 2016 Tomonori Kouya, All rights reserved */
/* Version 0.0: 2016-11-17 (Thu) First published          */
/**********************************************************/
#ifndef __TEMPLATE_TK_LINEAR_MP__
#define __TEMPLATE_TK_LINEAR_MP__

// return maximum value
#define MAX(i, j) (((i) > (j)) ? (i) : (j))

// 0-based index for square matrix
#define ZERO_INDEX(i, j, dim) ((i) * (dim) + (j))

// delete array (normal "delete" is not valid!)
#ifdef __MPREAL_H__
template <typename T = mpfr::mpreal> void delete_array(mpfr::mpreal array[], int dim)
{
	for(int i = 0; i < dim; i++)
		mpfr_clear(array[i].mpfr_ptr());

	delete [] array;
}
#endif // __MPREAL_H__

// delete array for double prec.
//template <typename T = double, dd_real, qd_real> void delete_array(double array[], int dim)
template <typename T> void delete_array(T array[], int dim)
{
	delete [] array;
}

// mymv : ret := A * x
// A : row-major order
template <typename T> void mymv(T ret[], T A[], T x[], int dim)
{
	int i, j;

	for(i = 0; i < dim; i++)
	{
		ret[i] = (T)0;
		for(j = 0; j < dim; j++)
			ret[i] += A[ZERO_INDEX(i, j, dim)] * x[j];
	}
}

// myaxpy : ret := alpha * x + y
template <typename T> void myaxpy(T ret[], T alpha, T x[], T y[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		ret[i] = alpha * x[i] + y[i];
}

// mydotp : x^T * y
template <typename T> T mydotp(T x[], T y[], int dim)
{
	int i;
	T ret = (T)0;

	for(i = 0; i < dim; i++)
		ret += x[i] * y[i];

	return ret;
}

// mynorm2 : ||x||_2
template <typename T> T mynorm2(T x[], int dim)
{
	int i;
	T ret = (T)0;

	for(i = 0; i < dim; i++)
		ret += x[i] * x[i];

	ret = (T)sqrt(ret);

	return ret;
}

// mycopy : x := y
template <typename T> void mycopy(T x[], T y[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		x[i] = y[i];
}

// set0 : x := 0
template <typename T> void set0(T x[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		x[i] = (T)0;
}

// relerr: return relative error of val
template <typename T> T get_relerr(T val, T true_val)
{
	if(true_val != (T)0)
		return abs((val - true_val) / true_val);
	else
		return abs(val - true_val);
}

// relerr_norm2: return relative error of val with norm2
template <typename T> T get_relerr_norm2(T vec[], T true_vec[], int dim)
{
	T true_vec_norm2, diff_norm2;
	T *diff;

	diff = new T[dim];

	// || true_vec ||_2
	true_vec_norm2 = mynorm2<T>(true_vec, dim);

	// diff := vec - true_vec
	myaxpy<T>(diff, (T)(-1), vec, true_vec, dim);
	diff_norm2 = mynorm2<T>(diff, dim);

	delete_array<T>(diff, dim);

	if(true_vec_norm2 != (T)0)
		return diff_norm2 / true_vec_norm2;
	else
		return diff_norm2;

}

// Solve A * x = b
template <typename T> int conjugate_gradient(T x[], T A[], T b[], int dim, T rel_tol, T abs_tol, int maxitimes)
{
	int itimes, i;

	// vectors
	T *r, *p, *w;

	// constant
	T alpha, beta, tmp_val, init_r_norm2, r_new_norm2;

	// Initialize
	r = new T[dim];
	p = new T[dim];
	w = new T[dim];

	// x_0 := 0
	set0<T>(x, dim);

	// r := b - A * x_0;
	mycopy<T>(r, b, dim); // x_0 = 0

	// init_r_norm2 = ||r||_2
	init_r_norm2 = mynorm2<T>(r, dim);

	// p := r
	mycopy<T>(p, r, dim);

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// w := A * p
		mymv<T>(w, A, p, dim);

		// alpha := (r, p) / (p, A * p)
		alpha = mydotp<T>(r, p, dim);
		tmp_val = mydotp<T>(p, w, dim);

		if(tmp_val == (T)0)
		{
			itimes = -1;
			break;
		}
		alpha /= tmp_val;

		// x := x + alpha * p
		myaxpy<T>(x, alpha, p, x, dim);

		// r_new := r - alpha * A * p
		beta = mydotp<T>(r, r, dim);
		myaxpy<T>(r, -alpha, w, r, dim);

		// beta := ||r_new||_2^2 / ||r||_2^2
		r_new_norm2 = mynorm2<T>(r, dim);
		beta = r_new_norm2 * r_new_norm2 / beta;

		// check residual
		if(r_new_norm2 <= rel_tol * init_r_norm2 + abs_tol)
			break;

		cout << setw(3) << itimes << " " << scientific << setprecision(3) << r_new_norm2 / init_r_norm2 << endl;

		// p := r + beta * p
		myaxpy<T>(p, beta, p, r, dim);

	}

	// clean
	delete_array<T>(r, dim);
	delete_array<T>(p, dim);
	delete_array<T>(w, dim);

	return itimes;
}

// LU decomposion of A
template <typename T> int LU(T A[], int dim, int pivot[])
{
	int i, j, k, max_j, tmp_index;
	T absmax_aji, abs_aji, pivot_aii;

	// initialize pivot vector
	for(i = 0; i < dim; i++)
		pivot[i] = i;

	// A decomposition
	for(i = 0; i < dim; i++)
	{
		// partial pivoting
		absmax_aji = abs(A[ZERO_INDEX(pivot[i], i, dim)]);
		max_j = i;
		for(j = i + 1; j < dim; j++)
		{
			abs_aji = A[ZERO_INDEX(pivot[j], i, dim)];
			if(absmax_aji < abs_aji)
			{
				max_j = j;
				absmax_aji = abs_aji;
			}
		}
		if(max_j != i)
		{
			tmp_index = pivot[max_j];
			pivot[max_j] = pivot[i];
			pivot[i] = tmp_index;
		}

		// select pivoted column
		pivot_aii = A[ZERO_INDEX(pivot[i], i, dim)];

		// error
		if(abs(pivot_aii) <= (T)0)
			return -1;

		for(j = i + 1; j < dim; j++)
		{
			A[ZERO_INDEX(pivot[j], i, dim)] /= pivot_aii;

			for(k = i + 1; k < dim; k++)
				A[ZERO_INDEX(pivot[j], k, dim)] -= A[ZERO_INDEX(pivot[j], i, dim)] * A[ZERO_INDEX(pivot[i], k, dim)];
		}
	}

	return 0;

}

// rowwise only
// solve LU * x = b in x -> b := x
template <typename T> int solve_LU_linear_eq(T x[], T LU[], T b[], int dim, int pivot[])
{
	int i, j;

	// x := b
	for(i = 0; i < dim; i++)
		x[i] = b[pivot[i]];

	// forward substitution
	for(j = 0; j < dim; j++)
	{
		for(i = j + 1; i < dim; i++)
			x[i] -= LU[ZERO_INDEX(pivot[i], j, dim)] * x[j];
	}

	// backward substitution
	for(i = dim - 1; i >= 0; i--)
	{
		for(j = i + 1; j < dim; j++)
			x[i] -= LU[ZERO_INDEX(pivot[i], j, dim)] * x[j];

		x[i] /= LU[ZERO_INDEX(pivot[i], i, dim)];
	}

	return 0;
}

// set matrix, true_x and b
template <typename T> void set_test_linear_eq(T A[], T true_x[], T b[], int dim)
{
	int i, j;

	// set A
	for(i = 0; i < dim; i++)
	{
		true_x[i] = (T)(i + 1);

		// Frank matrix
		//for(j = 0; j < dim; j++)
		//	A[ZERO_INDEX(i, j, dim)] = (T)((T)dim - (T)MAX(i, j));

		// Hilbert matrix
		for(j = 0; j < dim; j++)
			A[ZERO_INDEX(i, j, dim)] = (T)1 / (T)(i + j + 1);

	}

	// b := A * true_x
	mymv<T>(b, A, true_x, dim);

	// print the linear equation
/*	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
			cout << scientific << setprecision(32) << A[ZERO_INDEX(i, j, dim)] << " ";
		
		cout << "  x   " << scientific << setprecision(32) << b[i] << endl;
	}
*/
}


#endif // __TEMPLATE_TK_LINEAR_MP__
