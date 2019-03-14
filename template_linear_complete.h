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
/*void subst(mpfr::mpreal ret, double src)
{
		mpfr_set_d(ret.mpfr_ptr(), src, mpfr::mpreal::get_default_rnd());
}
*/
// double := mpfr
/*void subst(double *ret, mpfr::mpreal src)
{
		*ret = mpfr_get_d(src.mpfr_ptr(), mpfr::mpreal::get_default_rnd());
}
*/
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
			std::cerr << "Error: LU decomposition cannot be proceeded at " << i << "th pivot = " << pivot_aii << " !" << std::endl;
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
	int itimes, i;

	// vectors
	T *y, *zero_vec;
	T old_eig, tmp_val;

	// Initialize
	y = new T[dim];
	zero_vec = new T[dim];

	set0<T>(zero_vec, dim);

	// x_0 := 1
	set1<T>(eig_v, dim);

	// old_eig = 0;
	old_eig = (T)0;

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// y := A * x_old
		mymv<T>(y, A, eig_v, dim);

		// gamma := (y, x_old) / (x_old, x_old)
		eig_max = mydotp<T>(y, eig_v, dim);
		tmp_val = mydotp<T>(eig_v, eig_v, dim);

		if(tmp_val == (T)0)
		{
			itimes = -1;
			break;
		}
		eig_max /= tmp_val;

		// check convergence
		if(abs(eig_max - old_eig) <= rel_tol * abs(old_eig) + abs_tol)
			break;

		old_eig = eig_max;

		std::cout << std::setw(3) << itimes << " " << std::scientific << std::setprecision(10) << eig_max << std::endl;

		// eig_v := y / ||y||_2
		myaxpy<T>(eig_v, (T)1 / mynorm2<T>(y, dim), y, zero_vec, dim);

	}

	// clean
	delete_array<T>(y, dim);
	delete_array<T>(zero_vec, dim);

	return itimes;
}

// Inverse Power Method on R with CG
template <typename T> int inverse_power_method_cg(T &eig_min, T eig_v[], T A[], int dim, T rel_tol, T abs_tol, int maxitimes)
{
	int itimes, i;

	// vectors
	T *y, *zero_vec;
	T old_eig, tmp_val;

	// Initialize
	y = new T[dim];
	zero_vec = new T[dim];

	set0<T>(zero_vec, dim);

	// x_0 := 1
	set1<T>(eig_v, dim);

	// old_eig = 0;
	old_eig = (T)0;

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// y := A^(-1) * x_old
		conjugate_gradient<T>(y, A, eig_v, dim, rel_tol, abs_tol, maxitimes);

		// gamma := (x_old, x_old) / (y, x_old)
		eig_min = mydotp<T>(y, eig_v, dim);
		tmp_val = mydotp<T>(eig_v, eig_v, dim);

		if(eig_min == (T)0)
		{
			itimes = -1;
			break;
		}
		eig_min = tmp_val / eig_min;

		// check convergence
		if(abs(eig_min - old_eig) <= rel_tol * abs(old_eig) + abs_tol)
			break;

		old_eig = eig_min;

		std::cout << std::setw(3) << itimes << " " << std::scientific << std::setprecision(10) << eig_min << std::endl;

		// eig_v := y / ||y||_2
		myaxpy<T>(eig_v, (T)1 / mynorm2<T>(y, dim), y, zero_vec, dim);

	}

	// clean
	delete_array<T>(y, dim);
	delete_array<T>(zero_vec, dim);

	return itimes;
}

// Inverse Power Method on R with LU decomposition
template <typename T> int inverse_power_method_lu(T &eig_min, T eig_v[], T A[], int dim, T rel_tol, T abs_tol, int maxitimes)
{
	int itimes, i, *pivot;

	// vectors
	T *y, *zero_vec;
	T old_eig, tmp_val;

	// Initialize
	y = new T[dim];
	zero_vec = new T[dim];
	pivot = new int[dim];

	set0<T>(zero_vec, dim);

	// x_0 := 1
	set1<T>(eig_v, dim);

	// old_eig = 0;
	old_eig = (T)0;

	// LU decomposition
	LU<T>(A, dim, pivot);

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// y := A^(-1) * x_old
		solve_LU_linear_eq<T>(y, A, eig_v, dim, pivot);

		// gamma := (x_old, x_old) / (y, x_old)
		eig_min = mydotp<T>(y, eig_v, dim);
		tmp_val = mydotp<T>(eig_v, eig_v, dim);

		if(eig_min == (T)0)
		{
			itimes = -1;
			break;
		}
		eig_min = tmp_val / eig_min;

		// check convergence
		if(abs(eig_min - old_eig) <= rel_tol * abs(old_eig) + abs_tol)
			break;

		old_eig = eig_min;

		std::cout << std::setw(3) << itimes << " " << std::scientific << std::setprecision(10) << eig_min << std::endl;

		// eig_v := y / ||y||_2
		myaxpy<T>(eig_v, (T)1 / mynorm2<T>(y, dim), y, zero_vec, dim);

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
			std::cout << std::scientific << std::setprecision(32) << A[ZERO_INDEX(i, j, dim)] << " ";
		
		std::cout << "  x   " << std::scientific << std::setprecision(32) << b[i] << std::endl;
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
			std::cout << std::scientific << std::setprecision(16) << A[ZERO_INDEX(i, j, dim)] << " ";
		
		std::cout << "  x   " << std::scientific << std::setprecision(16) << b[i] << std::endl;
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

/* Referred from "A Numerical Comparison of Several Codition Number Estimators */
/* written by T.Matsuo, S.Masaaki, M.Mori (JSIAM vol.7, no.3 (1997) pp.307-319 */

/* coef_f_itr = rho_f(n) * cond_1(A) * eps_d / (1 - psi_f(n) * cond_1(A) * eps_s) */
template <typename T> T coef_f_itr(T cond1, T eps_d, T eps_s, int dim)
{
	T rho_f, psi_f;

	rho_f = sqrt((T)dim);
	psi_f = sqrt((T)dim);

	return rho_f * cond1 * eps_d / ((T)1 - psi_f * cond1 * eps_s);
}

/* Condition Number Estimator : cond_1(A) = ||A||_1 * ||A^(-1)||_1 */
/* Based on LINPACK argorithm */
/* norm1_orgmat: ||A||_1 */
/* lumat: LU(=A) decomposed matrix */
/* row_ch, col_ch: indeces for row and column */
template <typename T> T condest_dmatrix(T norm1_orgmat, T lumat[], T org_mat[], int dim, int pivot[])
{
	T condition_num, t, evec_i;
	T *zvec, *yvec, *xvec;
	long int i, j;
	long int row_i, col_i, *row_ch;

	/* Initialize */
	zvec = new T[dim];
	yvec = new T[dim];
	xvec = new T[dim];
	row_ch = new long int[dim];

	/* solve U^T * z = e */
	evec_i = (T)1;
	zvec[0] = evec_i / lumat[ZERO_INDEX(pivot[0], 0, dim)];
	for(i = 1; i < dim; i++)
	{
		row_i = pivot[i];
		col_i = i;

		/* t = sum^(i-1)_{j=0} u_ji * z_j */
		t = (T)0;
		for(j = 0; j <= (i - 1); j++)
			t += lumat[ZERO_INDEX(row_ch[j], col_i, dim)] * zvec[j];

		/* e[i] = sign(t) */
		if(t >= 0)
			evec_i = (T)1;
		else
			evec_i = (T)(-1);
		
		/* z[i] = (e[i] - t) / u[i][i] */
		zvec[i] = (evec_i - t) / lumat[ZERO_INDEX(row_i, col_i, dim)];
	}

	std::cout << "||z_vec||_1 = " << mynorm1(zvec) << std::endl;

	/* solve L^T (P * y) = z */
	/* L[i][i] = 1 */
	for(i = (dim - 1); i >= 0; i--)
	{
		row_i = row_ch[i];
		col_i = i;

		t = zvec[i];
		for(j = i + 1; j < dim; j++)
			t -= lumat[ZERO_INDEX(row_ch[j], col_i, dim)] * zvec[j];

		yvec[col_i] = t;
	}

	std::cout << "||y_vec||_1 = " << mynorm1(yvec, dim);

	/* solve A * x = y */
	solve_LU_linear_eq(xvec, lumat, yvec, dim, pivot);

	std::cout << "||x_vec||_1 = " << mynorm1(xvec, dim) << std::endl;

	/* cond_1(A) = ||A||_1 * ||A^(-1)||_1 \approx ||A||_1 * ||x||_1 / ||y||_1 */
	condition_num = mynorm1_mat(org_mat, dim) * mynorm1(xvec, dim) / mynorm1(yvec, dim);

	/* free */
	delete_array<T>(xvec, dim);
	delete_array<T>(yvec, dim);
	delete_array<T>(zvec, dim);
	delete row_ch;

	return condition_num;
}

// interative refinement with single & mpf_t mixed precision arithmetic
// solve x = a * b where known a in M_n(R) and b in R^n, unknown x in R^n
// unsigned long long_prec, short_prec; // long_prec > short_prec
template <typename L, typename S> int iterative_refinement(L x[], L a[], L b[], L rtol, L atol, int dim, int maxtimes)
{
	int i, itimes;
	int *af_ch;
	S *af;
	S *bf, *xf, *resf, *zf;
	L *res, *z;
	L tmp, norm_a, norm_x, norm_res, normalization_coef;

	// Initialize
	af = new S[dim * dim];
	bf = new S[dim];
	xf = new S[dim];

	res = new L[dim];

	resf = new S[dim];
	z = new L[dim];
	zf = new S[dim];
	af_ch = new int[dim];

	// norm_a := ||A||_F
	norm_a = mynorm2<L>(a, dim * dim);

	// Make short precision copy of A and b
	set_array(af, a, dim * dim);
	set_array(bf, b, dim);

	// Compute LU factorization in short precision
	// LU decomposition with partial pivoting
	LU<S>(af, dim, af_ch);

	// Apply back-solve in short precision with short precision factors
	solve_LU_linear_eq<S>(xf, af, bf, dim, af_ch);

	// Promote te solution from short precision to long precision
	set_array(x, xf, dim);

	// repeat iterative refinement process
	for(itimes = 0; itimes < maxtimes; itimes++)
	{
		// Compute residual in long precision
		// res = b - a * x
		mymv<L>(z, a, x, dim);
		myaxpy<L>(res, (L)(-1), z, b, dim);

		// until ||r_i||_2 < sqrt(n) * reps * ||A||_F * ||x_i||_2
		norm_x = mynorm2<L>(x, dim);
		norm_res = mynorm2<L>(res, dim);

		if(norm_res < sqrt((L)dim) * rtol * norm_a * norm_x + atol)
			break;

		// normalization: res := coef * res
		normalization_coef = (L)1 / norm_res;
		myscal<L>(res, normalization_coef, res, dim);

		// Demote the residual from long precision to short precison
		set_array(resf, res, dim);

		// Back-solve on short precision residual and short precision factors
		solve_LU_linear_eq<S>(zf, af, resf, dim, af_ch);

		// Promote the correction from short precision to long precision
		set_array(z, zf, dim);

		// reverse normalizationi
		myscal<L>(z, norm_res, z, dim);

		// Update solution in long precision
		myaxpy<L>(x, (L)1, x, z, dim);

		// for debug
		std::cout << itimes << ", " << std::setprecision(10) << norm_res << std::endl;
	}

	// if fail, retry in mpf_t precision
	if(itimes >= maxtimes)
	{
		// mpf_t precision
		LU<L>(a, dim, af_ch);
		//SolveMPFLS(x, a, b);
		solve_LU_linear_eq<L>(x, a, b, dim, af_ch);
	}

	// Clear
	delete_array<S>(af, dim);
	delete_array<S>(bf, dim);
	delete_array<S>(xf, dim);
	delete_array<S>(resf, dim);
	delete_array<L>(res, dim);
	delete_array<S>(zf, dim);
	delete_array<L>(z, dim);
	delete[] af_ch;

	return itimes;
}


#endif // __TEMPLATE_TK_LINEAR_MP__
