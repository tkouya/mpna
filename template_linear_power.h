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

// mpfr := double
void subst(mpfr::mpreal ret, double src)
{
		mpfr_set_d(ret.mpfr_ptr(), src, mpfr::mpreal::get_default_rnd());
}
// double := mpfr
void subst(double *ret, mpfr::mpreal src)
{
		*ret = mpfr_get_d(src.mpfr_ptr(), mpfr::mpreal::get_default_rnd());
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

// myscal : ret := alpha * x
template <typename T> void myscal(T ret[], T alpha, T x[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		ret[i] = alpha * x[i];
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

// mynorm1 : ||x||_1
template <typename T> T mynorm1(T x[], int dim)
{
	int i;
	T ret = (T)0;

	for(i = 0; i < dim; i++)
		ret += abs(x[i]);

	return ret;
}

// mynorm1 : ||x||_inf
template <typename T> T mynormi(T x[], int dim)
{
	int i;
	T ret = (T)abs(x[0]), abs_xi;

	for(i = 1; i < dim; i++)
	{
		abs_xi = abs(x[i]);
		if(ret < abs_xi)
			ret = abs_xi;
	}

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

// matrix norm1
template <typename T> T mynorm1_mat(T A[], int dim)
{
	int i, j;
	T ret, sum_abs_aij;

	ret = (T)0;
	for(j = 0; j < dim; j++)
	{
		sum_abs_aij = (T)0;
		for(i = 0; i < dim; i++)
			sum_abs_aij += abs(A[ZERO_INDEX(i, j, dim)]);

		if(ret < sum_abs_aij)
			ret = sum_abs_aij;
	}
}

// matrix norm_inf
template <typename T> T mynormi_mat(T A[], int dim)
{
	int i, j;
	T ret, sum_abs_aij;

	ret = (T)0;
	for(i = 0; i < dim; i++)
	{
		sum_abs_aij = (T)0;
		for(j = 0; j < dim; i++)
			sum_abs_aij += abs(A[ZERO_INDEX(i, j, dim)]);

		if(ret < sum_abs_aij)
			ret = sum_abs_aij;
	}
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

		//cout << setw(3) << itimes << " " << scientific << setprecision(3) << r_new_norm2 / init_r_norm2 << endl;

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
		absmax_aji = fabs(A[ZERO_INDEX(pivot[i], i, dim)]);
		max_j = i;
		for(j = i + 1; j < dim; j++)
		{
			abs_aji = fabs(A[ZERO_INDEX(pivot[j], i, dim)]);
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
		if(fabs(pivot_aii) <= (T)0)
		{
			std::cerr << "Error: LU decomposition cannot be proceeded at " << i << "th pivot = " << pivot_aii << " !" << endl;
			return -1;
		}

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


// set1 : x := 1
template <typename T> void set1(T x[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		x[i] = (T)1;
}

// Power Method on R
template <typename T> int power_method(T &eig_max, T eig_v[], T A[], int dim, T rel_tol, T abs_tol, int maxitimes)
{
	int itimes, i, nz_index;

	// vectors
	T *y, *zero_vec;
	T old_eig, norm_y;

	// Initialize
	y = new T[dim];
	zero_vec = new T[dim];

	set0<T>(zero_vec, dim);

	// eig_v_i := 1 -> eig_v / ||eig_v||_2
	set1<T>(eig_v, dim);
	norm_y = mynorm2<T>(eig_v, dim);
	myaxpy<T>(eig_v, (T)1 / norm_y, eig_v, zero_vec, dim);

	// old_eig = 0;
	old_eig = (T)0;

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// y := A * x_old
		mymv<T>(y, A, eig_v, dim);
		norm_y = mynorm2<T>(y, dim);

		// search nonzero element
		nz_index = -1;
		for(i = 0; i < dim; i++)
		{
			if(eig_v[i] != (T)0)
			{
				nz_index = i;
				break;
			}
		}
		if(nz_index == -1)
		{
			cerr << "ERROR in power_method: eigen vector is zero!" << endl;
			break;
		}

		// element-wise check
		eig_max = y[nz_index] / eig_v[nz_index];

		// check convergence
		if(abs(eig_max - old_eig) <= rel_tol * abs(old_eig) + abs_tol)
			break;

		old_eig = eig_max;

		cout << setw(3) << itimes << " " << scientific << setprecision(10) << eig_max << endl;

		// eig_v := y / ||y||_2
		myaxpy<T>(eig_v, (T)1 / norm_y, y, zero_vec, dim);

	}

	// clean
	delete_array<T>(y, dim);
	delete_array<T>(zero_vec, dim);

	return itimes;
}

// Inverse Power Method on R with CG
template <typename T> int inverse_power_method_cg(T &eig_min, T eig_v[], T A[], int dim, T rel_tol, T abs_tol, int maxitimes)
{
	int itimes, i, nz_index;

	// vectors
	T *y, *zero_vec;
	T old_eig, norm_y;

	// Initialize
	y = new T[dim];
	zero_vec = new T[dim];

	set0<T>(zero_vec, dim);

	// eig_v_i := 1 -> eig_v / ||eig_v||_2
	set1<T>(eig_v, dim);
	norm_y = mynorm2<T>(eig_v, dim);
	myaxpy<T>(eig_v, (T)1 / norm_y, eig_v, zero_vec, dim);

	// old_eig = 0;
	old_eig = (T)0;

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// y := A^(-1) * x_old
		conjugate_gradient<T>(y, A, eig_v, dim, rel_tol, abs_tol, maxitimes);

		norm_y = mynorm2<T>(y, dim);

		// search nonzero element
		nz_index = -1;
		for(i = 0; i < dim; i++)
		{
			if(y[i] != (T)0)
			{
				nz_index = i;
				break;
			}
		}
		if(nz_index == -1)
		{
			cerr << "ERROR in power_method: eigen vector is zero!" << endl;
			break;
		}

		// element-wise check
		eig_min = eig_v[nz_index] / y[nz_index];

		// check convergence
		if(abs(eig_min - old_eig) <= rel_tol * abs(old_eig) + abs_tol)
			break;

		old_eig = eig_min;

		cout << setw(3) << itimes << " " << scientific << setprecision(10) << eig_min << endl;

		// eig_v := y / ||y||_2
		myaxpy<T>(eig_v, (T)1 / norm_y, y, zero_vec, dim);

	}

	// clean
	delete_array<T>(y, dim);
	delete_array<T>(zero_vec, dim);

	return itimes;
}

// Inverse Power Method on R with LU decomposition
template <typename T> int inverse_power_method_lu(T &eig_min, T eig_v[], T A[], int dim, T rel_tol, T abs_tol, int maxitimes)
{
	int itimes, i, *pivot, nz_index;

	// vectors
	T *y, *zero_vec;
	T old_eig, norm_y;

	// Initialize
	y = new T[dim];
	zero_vec = new T[dim];
	pivot = new int[dim];

	set0<T>(zero_vec, dim);

	// eig_v_i := 1 -> eig_v / ||eig_v||_2
	set1<T>(eig_v, dim);
	norm_y = mynorm2<T>(eig_v, dim);
	myaxpy<T>(eig_v, (T)1 / norm_y, eig_v, zero_vec, dim);

	// old_eig = 0;
	old_eig = (T)0;

	// LU decomposition
	LU<T>(A, dim, pivot);

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// y := A^(-1) * x_old
		solve_LU_linear_eq<T>(y, A, eig_v, dim, pivot);

		norm_y = mynorm2<T>(y, dim);

		// search nonzero element
		nz_index = -1;
		for(i = 0; i < dim; i++)
		{
			if(y[i] != (T)0)
			{
				nz_index = i;
				break;
			}
		}
		if(nz_index == -1)
		{
			cerr << "ERROR in power_method: eigen vector is zero!" << endl;
			break;
		}

		// element-wise check
		eig_min = eig_v[nz_index] / y[nz_index];

		// check convergence
		if(abs(eig_min - old_eig) <= rel_tol * abs(old_eig) + abs_tol)
			break;

		old_eig = eig_min;

		cout << setw(3) << itimes << " " << scientific << setprecision(10) << eig_min << endl;

		// eig_v := y / ||y||_2
		myaxpy<T>(eig_v, (T)1 / norm_y, y, zero_vec, dim);

	}

	// clean
	delete_array<T>(y, dim);
	delete_array<T>(zero_vec, dim);
	delete pivot;

	return itimes;
}



// set matrix, true_x and b
template <typename T> void set_test_matrix(T A[], int dim)
{
	int i, j;

	// set A
	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
			A[ZERO_INDEX(i, j, dim)] = (T)((T)dim - (T)MAX(i, j));
	}

}

// set matrix, true_x and b
template <typename T> void set_test_linear_eq(T A[], T true_x[], T b[], int dim)
{
	int i, j;

	// set A
	for(i = 0; i < dim; i++)
	{
		true_x[i] = (T)(i + 1);
		for(j = 0; j < dim; j++)
			A[ZERO_INDEX(i, j, dim)] = (T)((T)dim - (T)MAX(i, j));
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

// set matrix, true_x and b
template <typename T> void set_test_linear_eq_hilbert(T A[], T true_x[], T b[], int dim)
{
	int i, j;

	// set A
	if(A != NULL)
	{
		for(i = 0; i < dim; i++)
		{
			for(j = 0; j < dim; j++)
				A[ZERO_INDEX(i, j, dim)] = (T)1/(T)(i + j + 1);
		}
	}

	// set true_x
	if(true_x != NULL)
	{
		for(i = 0; i < dim; i++)
			true_x[i] = (T)(i + 1);
	}

	// b := A * true_x
	if((b != NULL) && (A != NULL) && (true_x != NULL))
		mymv<T>(b, A, true_x, dim);

	// print the linear equation
/*	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
			cout << scientific << setprecision(16) << A[ZERO_INDEX(i, j, dim)] << " ";
		
		cout << "  x   " << scientific << setprecision(16) << b[i] << endl;
	}
*/
}


// substitute the same type
template <typename T> void set_array(T ret_array[], T org_array[], int dim)
{
	for(int i = 0; i < dim; i++)
		ret_array[i] = org_array[i];
}

#ifdef _QD_QD_REAL_H

// conversion of data types
// double -> dd_real
void set_array(dd_real ret_array[], double org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
	{
		ret_array[i].x[0] = org_array[i];
		ret_array[i].x[1] = (double)0.0;
	}
}

// conversion of data types
// double -> qd_real
void set_array(qd_real ret_array[], double org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
	{
		ret_array[i].x[0] = org_array[i];
		ret_array[i].x[1] = (double)0.0;
		ret_array[i].x[2] = (double)0.0;
		ret_array[i].x[3] = (double)0.0;
	}
}

// conversion of data types
// dd_real -> double
void set_array(double ret_array[], dd_real org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
		ret_array[i] = org_array[i].x[0];

}

// conversion of data types
// qd_real -> double
void set_array(double ret_array[], qd_real org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
		ret_array[i] = org_array[i].x[0];

}

// conversion of data types
// dd_real -> qd_real
void set_array(qd_real ret_array[], dd_real org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
	{
		ret_array[i].x[0] = org_array[i].x[0];
		ret_array[i].x[1] = org_array[i].x[1];
		ret_array[i].x[2] = 0.0;
		ret_array[i].x[3] = 0.0;
	}
}

// conversion of data types
// qd_real -> dd_real
void set_array(dd_real ret_array[], qd_real org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
	{
		ret_array[i].x[0] = org_array[i].x[0];
		ret_array[i].x[1] = org_array[i].x[1];
	}
}

#endif // _QD_QD_REAL_H

#ifdef __MPREAL_H__
// conversion of data types
// mpreal -> double
void set_array(double ret_array[], mpfr::mpreal org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
		ret_array[i] = mpfr_get_d(org_array[i].mpfr_ptr(), MPFR_RNDN);

}

// conversion of data types
// double -> mpreal
void set_array(mpfr::mpreal ret_array[], double org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
		ret_array[i] = (mpfr::mpreal)org_array[i];

}

#ifdef _QD_QD_REAL_H

// mpfr_[get, set]_[dd, qd]
#include "mpfr_dd_qd.h"

// conversion of data types
// mpreal -> dd_real
void set_array(dd_real ret_array[], mpfr::mpreal org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
		mpfr_get_dd(ret_array[i].x, org_array[i].mpfr_ptr(), MPFR_RNDN);

}

// conversion of data types
// mpreal -> qd_real
void set_array(qd_real ret_array[], mpfr::mpreal org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
		mpfr_get_qd(ret_array[i].x, org_array[i].mpfr_ptr(), MPFR_RNDN);

}

// conversion of data types
// dd_real -> mpreal
void set_array(mpfr::mpreal ret_array[], dd_real org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
		mpfr_set_dd(ret_array[i].mpfr_ptr(), org_array[i].x, MPFR_RNDN);

}

// conversion of data types
// qd_real -> mpreal
void set_array(mpfr::mpreal ret_array[], qd_real org_array[], int dim)
{
	// double -> dd_real, qd_real, mpreal
	for(int i = 0; i < dim; i++)
		mpfr_set_qd(ret_array[i].mpfr_ptr(), org_array[i].x, MPFR_RNDN);

}

#endif // _QD_QD_REAL_H

#endif // __MPREAL_H__


#endif // __TEMPLATE_TK_LINEAR_MP__
