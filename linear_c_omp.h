//******************************************************************************
// linear_c_omp.h : C functions for Multiple precision Parallelized 
//                                                linear computation with OpenMP
//                              S                  upporting double and MPFR/GMP
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
#ifndef __LINEAR_TK_C_MP_OMP__
#define __LINEAR_TK_C_MP_OMP__

#include <omp.h> // OpenMP header file

// return maximum value
#define MAX(i, j) (((i) > (j)) ? (i) : (j))

// 0-based index for square matrix
#define ZERO_INDEX(i, j, dim) ((i) * (dim) + (j))

///////////////////////////////////
// double precision
///////////////////////////////////

// mymv : ret := A * x
// A : row-major order
void d_mymv(double ret[], double A[], double x[], int dim)
{
	int i, j;

	#pragma omp parallel for private(i, j) shared(ret, A, x, dim)
	for(i = 0; i < dim; i++)
	{
		ret[i] = 0.0;

		//#pragma omp parallel for private(j) shared(ret, A, x, dim)
		for(j = 0; j < dim; j++)
			ret[i] += A[ZERO_INDEX(i, j, dim)] * x[j];
	}
}

// myaxpy : ret := alpha * x + y
void d_myaxpy(double ret[], double alpha, double x[], double y[], int dim)
{
	int i;

	#pragma omp parallel for private(i) shared(dim, ret, alpha, x, y)
	for(i = 0; i < dim; i++)
		ret[i] = alpha * x[i] + y[i];
}

// mydotp : x^T * y
double d_mydotp(double x[], double y[], int dim)
{
	int i, thread_index;
	double ret = 0.0;

	#pragma omp parallel for private(i) shared(x, y, dim) reduction(+:ret)
	for(i = 0; i < dim; i++)
		ret += x[i] * y[i];

	return ret;
}

// mynorm2 : ||x||_2
double d_mynorm2(double x[], int dim)
{
	int i;
	double ret = 0.0;

	#pragma omp parallel for private(i) shared(x, dim) reduction(+:ret)
	for(i = 0; i < dim; i++)
		ret += x[i] * x[i];

	ret = sqrt(ret);

	return ret;
}

// mycopy : x := y
void d_mycopy(double x[], double y[], int dim)
{
	int i;

	#pragma omp parallel for private(i)
	for(i = 0; i < dim; i++)
		x[i] = y[i];
}

// set0 : x := 0
void d_set0(double x[], int dim)
{
	int i;

	#pragma omp parallel for private(i)
	for(i = 0; i < dim; i++)
		x[i] = 0.0;
}

// relerr: return relative error of val
double get_d_relerr(double val, double true_val)
{
	if(true_val != 0.0)
		return fabs((val - true_val) / true_val);
	else
		return fabs(val - true_val);
}

// relerr_norm2: return relative error of val with norm2
double get_d_relerr_norm2(double vec[], double true_vec[], int dim)
{
	double true_vec_norm2, diff_norm2;
	double *diff;

	diff = (double *)calloc(dim, sizeof(double));

	// || true_vec ||_2
	true_vec_norm2 = d_mynorm2(true_vec, dim);

	// diff := vec - true_vec
	d_myaxpy(diff, -1.0, vec, true_vec, dim);
	diff_norm2 = d_mynorm2(diff, dim);

	free(diff);

	if(true_vec_norm2 != 0.0)
		return diff_norm2 / true_vec_norm2;
	else
		return diff_norm2;

}

// Solve A * x = b
int d_conjugate_gradient(double x[], double A[], double b[], int dim, double rel_tol, double abs_tol, int maxitimes)
{
	int itimes, i;

	// vectors
	double *r, *p, *w;

	// constant
	double alpha, beta, tmp_val, init_r_norm2, r_new_norm2;

	// Initialize
	r = (double *)calloc(dim, sizeof(double));
	p = (double *)calloc(dim, sizeof(double));
	w = (double *)calloc(dim, sizeof(double));

	// x_0 := 0
	d_set0(x, dim);

	// r := b - A * x_0;
	d_mycopy(r, b, dim); // x_0 = 0

	// init_r_norm2 = ||r||_2
	init_r_norm2 = d_mynorm2(r, dim);

	// p := r
	d_mycopy(p, r, dim);

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// w := A * p
		d_mymv(w, A, p, dim);

		// alpha := (r, p) / (p, A * p)
		alpha = d_mydotp(r, p, dim);
		tmp_val = d_mydotp(p, w, dim);

		if(tmp_val == 0.0)
		{
			itimes = -1;
			break;
		}
		alpha /= tmp_val;

		// x := x + alpha * p
		d_myaxpy(x, alpha, p, x, dim);

		// r_new := r - alpha * A * p
		beta = d_mydotp(r, r, dim);
		d_myaxpy(r, -alpha, w, r, dim);

		// beta := ||r_new||_2^2 / ||r||_2^2
		r_new_norm2 = d_mynorm2(r, dim);
		beta = r_new_norm2 * r_new_norm2 / beta;

		// check residual
		if(r_new_norm2 <= rel_tol * init_r_norm2 + abs_tol)
			break;

//		printf("%3d %10.3e\n", itimes, r_new_norm2 / init_r_norm2);

		// p := r + beta * p
		d_myaxpy(p, beta, p, r, dim);

	}

	// clean
	free(r);
	free(p);
	free(w);

	return itimes;
}

// set matrix, true_x and b
void set_test_d_linear_eq(double A[], double true_x[], double b[], int dim)
{
	int i, j;

	// set A
	#pragma omp parallel for private(i, j) shared(true_x, A, dim)
	for(i = 0; i < dim; i++)
	{
		true_x[i] = (double)(i + 1);
		for(j = 0; j < dim; j++)
			A[ZERO_INDEX(i, j, dim)] = (double)(dim - MAX(i, j));
	}

	// b := A * true_x
	d_mymv(b, A, true_x, dim);

	// print the linear equation
/*	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
			printf("%25.17e ", A[ZERO_INDEX(i, j, dim)]);
		
		printf("  x   %25.17e\n", b[i]);
	}
*/
}

// LU decomposion of A
int d_LU(double A[], int dim, int pivot[])
{
	int i, j, k, max_j, tmp_index;
	double absmax_aji, abs_aji, pivot_aii;

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
		if(fabs(pivot_aii) <= (double)0)
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
int d_solve_LU_linear_eq(double x[], double LU[], double b[], int dim, int pivot[])
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

///////////////////////////////////
// MPFR/GMP
///////////////////////////////////
#ifdef __MPFR_H

// default rounding mode = round to nearest
mpfr_rnd_t _tk_default_rmode = MPFR_RNDN;

// DBL_MAX
#include <float.h>

// FALSE(=0), TRUE(=1)
//#define MPFR_IS_SINGULAR(x) (MPFR_EXP(x) <= MPFR_EXP_INF)
int mpfr_is_singular(mpfr_srcptr x)
{
	//return (mpfr_get_exp(x) <= MPFR_EXP_INF);
	return (mpfr_get_exp(x) <= mpfr_get_emin());
}

/* generic code */
// ret_dd[2] = ret_dd[high == 0], ret_dd[low == 1]
void mpfr_get_dd (double ret_dd[2], mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
	//if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
	if (mpfr_is_singular(x))
	{
		//return (long double) mpfr_get_d (x, rnd_mode);
		ret_dd[0] = mpfr_get_d(x, rnd_mode);
		ret_dd[1] = (double)0.0;
	}
  else /* now x is a normal non-zero number */
	{
      // long double r; /* result */
      double s; /* part of result */
      //MPFR_SAVE_EXPO_DECL (expo);

      //MPFR_SAVE_EXPO_MARK (expo);

//#if defined(HAVE_LDOUBLE_MAYBE_DOUBLE_DOUBLE)
      //if (MPFR_LDBL_MANT_DIG == 106)
      //{
          /* Assume double-double format (as found with the PowerPC ABI).
             The generic code below isn't used because numbers with
             precision > 106 would not be supported. */
          s = mpfr_get_d (x, MPFR_RNDN); /* high part of x */
          /* Let's first consider special cases separately. The test for
             infinity is really needed to avoid a NaN result. The test
             for NaN is mainly for optimization. The test for 0 is useful
             to get the correct sign (assuming mpfr_get_d supports signed
             zeros on the implementation). */
          //if (s == 0 || DOUBLE_ISNAN (s) || DOUBLE_ISINF (s))
          if (s == 0 || isnan(s) || isinf(s))
					{
            // r = (long double) s;
						ret_dd[0] = s;
						ret_dd[1] = 0;
					}
          else
					{
              mpfr_t y, z;

              mpfr_init2 (y, mpfr_get_prec (x));
              mpfr_init2 (z, 53); /* keep the precision small */
              mpfr_set_d (z, s, MPFR_RNDN);  /* exact */

							// y := x - z = x - s
              mpfr_sub (y, x, z, MPFR_RNDN); /* exact */

              /* Add the second part of y (in the correct rounding mode). */
              //r = (long double) s + (long double) mpfr_get_d (y, rnd_mode);
							ret_dd[0] = s;
							ret_dd[1] = mpfr_get_d(y, rnd_mode);

              mpfr_clear (z);
              mpfr_clear (y);
          }
      //}
      //MPFR_SAVE_EXPO_FREE (expo);
      //return r;
    }
}

/* double-double code */
int mpfr_set_dd (mpfr_ptr r, double d[2], mpfr_rnd_t rnd_mode)
{
  mpfr_t t, u;
  int inexact, shift_exp;
  double h, l;
//  MPFR_SAVE_EXPO_DECL (expo);

  /* Check for NAN */
  //ONGDOUBLE_NAN_ACTION (d, goto nan);
	if(isnan(d[0])) goto nan;

  /* Check for INF */
//  if (d > MPFR_LDBL_MAX)
  if (d[0] > DBL_MAX)
    {
      mpfr_set_inf (r, 1);
      return 0;
    }
  //else if (d < -MPFR_LDBL_MAX)
  else if (d[0] < -DBL_MAX)
    {
      mpfr_set_inf (r, -1);
      return 0;
    }
  /* Check for ZERO */
  else if (d[0] == 0.0)
    return mpfr_set_d (r, (double) d[0], rnd_mode);

 // if (d >= (long double) MPFR_LDBL_MAX || d <= (long double) -MPFR_LDBL_MAX)
 //   h = (d >= (long double) MPFR_LDBL_MAX) ? MPFR_LDBL_MAX : -MPFR_LDBL_MAX;
  if (d[0] >= DBL_MAX || d[0] <= -DBL_MAX)
    h = (d[0] >= DBL_MAX) ? DBL_MAX : -DBL_MAX;
  else
    h = (double) d[0]; /* should not overflow */
  //l = (double) (d - (long double) h);
	l = d[1];

  //MPFR_SAVE_EXPO_MARK (expo);

  mpfr_init2 (t, 53);
  mpfr_init2 (u, 53);

  inexact = mpfr_set_d (t, h, MPFR_RNDN);
  //MPFR_ASSERTN(inexact == 0);
  inexact = mpfr_set_d (u, l, MPFR_RNDN);
  //MPFR_ASSERTN(inexact == 0);
  inexact = mpfr_add (r, t, u, rnd_mode);

  mpfr_clear (t);
  mpfr_clear (u);

//  MPFR_SAVE_EXPO_FREE (expo);
  inexact = mpfr_check_range (r, inexact, rnd_mode);
  return inexact;

 nan:
  mpfr_set_nan(r); //MPFR_SET_NAN(r);
  //return MPFR_RET_NAN;
	return 0;
}

/* generic code */
// ret_qd[4] = ret_qd[high == 0], ret_qd[1], ret_qd[2], ret_qd[3]
void mpfr_get_qd (double ret_qd[4], mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
	int i;
	double s;

	//if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
	if (mpfr_is_singular(x))
	{
		//return (long double) mpfr_get_d (x, rnd_mode);
		ret_qd[0] = mpfr_get_d(x, rnd_mode);
		ret_qd[1] = (double)0.0;
		ret_qd[2] = (double)0.0;
		ret_qd[3] = (double)0.0;
	}
	else /* now x is a normal non-zero number */
	{
//#if defined(HAVE_LDOUBLE_MAYBE_DOUBLE_DOUBLE)
      //if (MPFR_LDBL_MANT_DIG == 106)
      //{
          /* Assume double-double format (as found with the PowerPC ABI).
             The generic code below isn't used because numbers with
             precision > 106 would not be supported. */
          s = mpfr_get_d (x, MPFR_RNDN); /* high part of x */
          /* Let's first consider special cases separately. The test for
             infinity is really needed to avoid a NaN result. The test
             for NaN is mainly for optimization. The test for 0 is useful
             to get the correct sign (assuming mpfr_get_d supports signed
             zeros on the implementation). */
          //if (s == 0 || DOUBLE_ISNAN (s) || DOUBLE_ISINF (s))
          if (s == 0 || isnan(s) || isinf(s))
					{
            // r = (long double) s;
						ret_qd[0] = s;
						ret_qd[1] = 0;
						ret_qd[1] = 0;
						ret_qd[1] = 0;
					}
          else
					{
              mpfr_t y, z[3];

              mpfr_init2 (y, mpfr_get_prec (x));
              mpfr_init2 (z[0], 53); /* keep the precision small */
              mpfr_init2 (z[1], 53); /* keep the precision small */
              mpfr_init2 (z[2], 53); /* keep the precision small */

							// y[0] := x - z = x - s
              mpfr_set_d (z[0], s, MPFR_RNDN);  /* exact */
              mpfr_sub (y, x, z[0], MPFR_RNDN); /* exact */
							ret_qd[0] = s;

							// y[1] := y[0] - z = x - s
 		          s = mpfr_get_d (y, MPFR_RNDN); /* high part of y */
	            mpfr_set_d (z[1], s, MPFR_RNDN);  /* exact */
              mpfr_sub (y, x, z[0], MPFR_RNDN); /* exact */
              mpfr_sub (y, y, z[1], MPFR_RNDN); /* exact */
							ret_qd[1] = s;

							// y[2] := y[1] - z = x - s
 		          s = mpfr_get_d (y, MPFR_RNDN); /* high part of y */
	            mpfr_set_d (z[2], s, MPFR_RNDN);  /* exact */
              mpfr_sub (y, x, z[0], MPFR_RNDN); /* exact */
              mpfr_sub (y, y, z[1], MPFR_RNDN); /* exact */
              mpfr_sub (y, y, z[2], MPFR_RNDN); /* exact */
							ret_qd[2] = s;

							// y[3] := x - z = x - s
							ret_qd[3] = mpfr_get_d(y, rnd_mode);

              mpfr_clear (z[0]);
              mpfr_clear (z[1]);
              mpfr_clear (z[2]);
              mpfr_clear (y);
          }
      //}
      //MPFR_SAVE_EXPO_FREE (expo);
      //return r;
    }
}

/* double-double code */
int mpfr_set_qd (mpfr_ptr r, double d[4], mpfr_rnd_t rnd_mode)
{
  mpfr_t t, u, v, w;
  int inexact, shift_exp;
  double h, l[3];
//  MPFR_SAVE_EXPO_DECL (expo);

  /* Check for NAN */
  //ONGDOUBLE_NAN_ACTION (d, goto nan);
	if(isnan(d[0])) goto nan;

  /* Check for INF */
//  if (d > MPFR_LDBL_MAX)
  if (d[0] > DBL_MAX)
    {
      mpfr_set_inf (r, 1);
      return 0;
    }
  //else if (d < -MPFR_LDBL_MAX)
  else if (d[0] < -DBL_MAX)
    {
      mpfr_set_inf (r, -1);
      return 0;
    }
  /* Check for ZERO */
  else if (d[0] == 0.0)
    return mpfr_set_d (r, (double) d[0], rnd_mode);

 // if (d >= (long double) MPFR_LDBL_MAX || d <= (long double) -MPFR_LDBL_MAX)
 //   h = (d >= (long double) MPFR_LDBL_MAX) ? MPFR_LDBL_MAX : -MPFR_LDBL_MAX;
  if (d[0] >= DBL_MAX || d[0] <= -DBL_MAX)
    h = (d[0] >= DBL_MAX) ? DBL_MAX : -DBL_MAX;
  else
    h = (double) d[0]; /* should not overflow */
  //l = (double) (d - (long double) h);
	l[0] = d[1];
	l[1] = d[2];
	l[2] = d[3];

  //MPFR_SAVE_EXPO_MARK (expo);

  mpfr_init2 (t, 53);
  mpfr_init2 (u, 53);
  mpfr_init2 (v, 53);
  mpfr_init2 (w, 53);

  inexact = mpfr_set_d (t, h, MPFR_RNDN);
  inexact = mpfr_set_d (u, l[0], MPFR_RNDN);
  inexact = mpfr_set_d (v, l[1], MPFR_RNDN);
  inexact = mpfr_set_d (w, l[2], MPFR_RNDN);

  inexact = mpfr_add (r, t, u, rnd_mode);
  inexact = mpfr_add (r, r, v, rnd_mode);
  inexact = mpfr_add (r, r, w, rnd_mode);

  mpfr_clear (t);
  mpfr_clear (u);
  mpfr_clear (v);
  mpfr_clear (w);

//  MPFR_SAVE_EXPO_FREE (expo);
  inexact = mpfr_check_range (r, inexact, rnd_mode);
  return inexact;

 nan:
  mpfr_set_nan(r); //MPFR_SET_NAN(r);
  //return MPFR_RET_NAN;
	return 0;
}

// mymv : ret := A * x
// A : row-major order
void mpfr_mymv(mpfr_t ret[], mpfr_t A[], mpfr_t x[], int dim)
{
	int i, j;

	#pragma omp parallel for private(i, j) shared(dim, ret, A, x, _tk_default_rmode)
	for(i = 0; i < dim; i++)
	{
		mpfr_set_ui(ret[i], 0UL, _tk_default_rmode);

		for(j = 0; j < dim; j++)
			mpfr_fma(ret[i], A[ZERO_INDEX(i, j, dim)], x[j], ret[i], _tk_default_rmode);
	}
}

// myaxpy : ret := alpha * x + y
void mpfr_myaxpy(mpfr_t ret[], mpfr_t alpha, mpfr_t x[], mpfr_t y[], int dim)
{
	int i;

	#pragma omp parallel for private(i) shared(dim, _tk_default_rmode)
	for(i = 0; i < dim; i++)
		mpfr_fma(ret[i], alpha, x[i], y[i], _tk_default_rmode);
}

// mydotp : x^T * y
void mpfr_mydotp(mpfr_t ret, mpfr_t x[], mpfr_t y[], int dim)
{
	unsigned long prec;
	mpfr_t tmp_val;
	int i;

	prec = mpfr_get_prec(ret);
	mpfr_set_ui(ret, 0UL, _tk_default_rmode);

	#pragma omp parallel private(tmp_val) shared(prec, dim, ret, _tk_default_rmode)
	{
		mpfr_init2(tmp_val, prec);
		mpfr_set_ui(tmp_val, 0UL, _tk_default_rmode);

		#pragma omp for private(i)
		for(i = 0; i < dim; i++)
			mpfr_fma(tmp_val, x[i], y[i], tmp_val, _tk_default_rmode);

		#pragma omp critical
		mpfr_add(ret, ret, tmp_val, _tk_default_rmode);

		mpfr_clear(tmp_val);

	} // omp parallel

	return;
}

// mynorm2 : ||x||_2
void mpfr_mynorm2(mpfr_t ret, mpfr_t x[], int dim)
{
	unsigned long prec;
	mpfr_t tmp_val;
	int i;

	prec = mpfr_get_prec(ret);
	mpfr_set_ui(ret, 0UL, _tk_default_rmode);

	#pragma omp parallel private(tmp_val) shared(prec, dim, ret, _tk_default_rmode)
	{
		mpfr_init2(tmp_val, prec);
		mpfr_set_ui(tmp_val, 0UL, _tk_default_rmode);

		#pragma omp for private(i)
		for(i = 0; i < dim; i++)
			mpfr_fma(tmp_val, x[i], x[i], tmp_val, _tk_default_rmode);

		#pragma omp critical
		mpfr_add(ret, ret, tmp_val, _tk_default_rmode);

		mpfr_clear(tmp_val);

	}

	mpfr_sqrt(ret, ret, _tk_default_rmode);

	return;
}

// mycopy : x := y
void mpfr_mycopy(mpfr_t x[], mpfr_t y[], int dim)
{
	int i;

	#pragma omp parallel for private(i) shared(dim, x, y, _tk_default_rmode)
	for(i = 0; i < dim; i++)
		mpfr_set(x[i], y[i], _tk_default_rmode);
}

// set0 : x := 0
void mpfr_set0(mpfr_t x[], int dim)
{
	int i;

	#pragma omp parallel for private(i) shared(dim, x, _tk_default_rmode)
	for(i = 0; i < dim; i++)
		mpfr_set_ui(x[i], 0UL, _tk_default_rmode);

}

// relerr: return relative error of val
void get_mpfr_relerr(mpfr_t ret, mpfr_t val, mpfr_t true_val)
{
	// ret := val - true_val
	mpfr_sub(ret, val, true_val, _tk_default_rmode);

	// ret := (val - true_val) / true_val
	if(mpfr_cmp_ui(true_val, 0UL) != 0)
		mpfr_div(ret, ret, true_val, _tk_default_rmode);

	// ret := |ret|
	mpfr_abs(ret, ret, _tk_default_rmode);	

	return;
}

// mpfr_init_array: initialize array by mpfr_init
void mpfr_init_array(mpfr_t array[], int dim)
{
	int i;

	#pragma omp parallel for private(i) shared(dim, array)
	for(i = 0; i < dim; i++)
		mpfr_init(array[i]);
}

// mpfr_init2_array: initialize array by mpfr_init2
void mpfr_init2_array(mpfr_t array[], int dim, unsigned long prec)
{
	int i;

	#pragma omp parallel for private(i) shared(dim, array)
	for(i = 0; i < dim; i++)
		mpfr_init2(array[i], prec);
}

// mpfr_clear_array: clear array by mpfr_clear
void mpfr_clear_array(mpfr_t array[], int dim)
{
	int i;

	#pragma omp parallel for private(i) shared(dim, array)
	for(i = 0; i < dim; i++)
		mpfr_init(array[i]);
}

// mpfr_get_max_prec_array: get max prec in array
unsigned long mpfr_get_max_prec_array(mpfr_t array[], int dim)
{
	int i;
	unsigned long max_prec, prec;

	max_prec = mpfr_get_prec(array[0]);

	#pragma omp parallel for private(i, prec) reduction(max: max_prec)
	for(i = 1; i < dim; i++)
	{
		prec = mpfr_get_prec(array[i]);
		if(prec > max_prec)
			max_prec = prec;
	}

	return max_prec;
}

// relerr_norm2: return relative error of val with norm2
void get_mpfr_relerr_norm2(mpfr_t ret, mpfr_t vec[], mpfr_t true_vec[], int dim)
{
	unsigned long prec;
	mpfr_t true_vec_norm2, diff_norm2, minus_one;
	mpfr_t *diff;

	// initialize
	prec = mpfr_get_prec(ret);
	mpfr_init2(true_vec_norm2, prec);
	mpfr_init2(diff_norm2, prec);
	mpfr_init2(minus_one, prec); mpfr_set_si(minus_one, -1, _tk_default_rmode);

	diff = (mpfr_t *)calloc(dim, sizeof(mpfr_t));
	mpfr_init2_array(diff, dim, prec);

	// || true_vec ||_2
	mpfr_mynorm2(true_vec_norm2, true_vec, dim);

	// diff := vec - true_vec
	mpfr_myaxpy(diff, minus_one, vec, true_vec, dim);
	mpfr_mynorm2(diff_norm2, diff, dim);

	if(mpfr_cmp_ui(true_vec_norm2, 0UL) != 0)
		mpfr_div(ret, diff_norm2, true_vec_norm2, _tk_default_rmode);
	else
		mpfr_set(ret, diff_norm2, _tk_default_rmode);

	// free
	mpfr_clear(true_vec_norm2);
	mpfr_clear(diff_norm2);
	mpfr_clear(minus_one);
	mpfr_clear_array(diff, dim);

	return;
}

// Solve A * x = b
int mpfr_conjugate_gradient(mpfr_t x[], mpfr_t A[], mpfr_t b[], int dim, mpfr_t rel_tol, mpfr_t abs_tol, int maxitimes)
{
	unsigned long prec;
	int itimes, i;

	// vectors
	mpfr_t *r, *p, *w;

	// constant
	mpfr_t alpha, beta, tmp_val, init_r_norm2, r_new_norm2;

	// Initialize
	prec = mpfr_get_max_prec_array(x, dim); // x's precision is locally used

	mpfr_inits2(prec, alpha, beta, tmp_val, init_r_norm2, r_new_norm2, (mpfr_ptr)NULL);

	r = (mpfr_t *)calloc(dim, sizeof(mpfr_t));
	p = (mpfr_t *)calloc(dim, sizeof(mpfr_t));
	w = (mpfr_t *)calloc(dim, sizeof(mpfr_t));

	mpfr_init2_array(r, dim, prec);
	mpfr_init2_array(p, dim, prec);
	mpfr_init2_array(w, dim, prec);

	// x_0 := 0
	mpfr_set0(x, dim);

	// r := b - A * x_0;
	mpfr_mycopy(r, b, dim); // x_0 = 0

	// init_r_norm2 = ||r||_2
	mpfr_mynorm2(init_r_norm2, r, dim);

	// p := r
	mpfr_mycopy(p, r, dim);

	// main loop
	for(itimes = 0; itimes < maxitimes; itimes++)
	{
		// w := A * p
		mpfr_mymv(w, A, p, dim);

		// alpha := (r, p) / (p, A * p)
		mpfr_mydotp(alpha, r, p, dim);
		mpfr_mydotp(tmp_val, p, w, dim);

		if(mpfr_cmp_ui(tmp_val, 0UL) == 0)
		{
			itimes = -1;
			break;
		}
		//alpha /= tmp_val;
		mpfr_div(alpha, alpha, tmp_val, _tk_default_rmode);

		// x := x + alpha * p
		mpfr_myaxpy(x, alpha, p, x, dim);

		// r_new := r - alpha * A * p
		mpfr_mydotp(beta, r, r, dim);
		mpfr_neg(alpha, alpha, _tk_default_rmode); // alpha := -alpha
		mpfr_myaxpy(r, alpha, w, r, dim);

		// beta := ||r_new||_2^2 / ||r||_2^2
		//beta = r_new_norm2 * r_new_norm2 / beta;
		mpfr_mynorm2(r_new_norm2, r, dim);
		mpfr_sqr(tmp_val, r_new_norm2, _tk_default_rmode);
		mpfr_div(beta, tmp_val, beta, _tk_default_rmode);

		// check residual
		//if(r_new_norm2 <= rel_tol * init_r_norm2 + abs_tol)
		mpfr_mul(tmp_val, rel_tol, init_r_norm2, _tk_default_rmode);
		mpfr_add(tmp_val, tmp_val, abs_tol, _tk_default_rmode);
		if(mpfr_cmp(r_new_norm2, tmp_val) <= 0)
			break;

		// tmp_val := r_new_norm2 / init_r_norm2
		mpfr_div(tmp_val, r_new_norm2, init_r_norm2, _tk_default_rmode);
		//mpfr_printf("%3d %10.3RNe\n", itimes, tmp_val);

		// p := r + beta * p
		mpfr_myaxpy(p, beta, p, r, dim);

	}

	// clean
	mpfr_clear_array(r, dim);
	mpfr_clear_array(p, dim);
	mpfr_clear_array(w, dim);
	free(r);
	free(p);
	free(w);

	mpfr_clears(alpha, beta, tmp_val, init_r_norm2, r_new_norm2, (mpfr_ptr)0);

	return itimes;
}

// set matrix, true_x and b
void set_test_mpfr_linear_eq(mpfr_t A[], mpfr_t true_x[], mpfr_t b[], int dim)
{
	int i, j;

	// set A
	#pragma omp parallel for private(i, j) shared(dim, true_x, A, _tk_default_rmode)
	for(i = 0; i < dim; i++)
	{
		mpfr_set_ui(true_x[i], (unsigned long)(i + 1), _tk_default_rmode);
		for(j = 0; j < dim; j++)
			mpfr_set_ui(A[ZERO_INDEX(i, j, dim)], (unsigned long)(dim - MAX(i, j)), _tk_default_rmode);
	}

	// b := A * true_x
	mpfr_mymv(b, A, true_x, dim);

	// print the linear equation
/*	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
			mpfr_printf("%25.17RNe ", A[ZERO_INDEX(i, j, dim)]);
		
		mpfr_printf("  x   %25.17RNe\n", b[i]);
	}
*/
}

// LU decomposion of A
int mpfr_LU(mpfr_t A[], int dim, int pivot[])
{
	int i, j, k, max_j, tmp_index;
	mpfr_t absmax_aji, abs_aji, pivot_aii, abs_pivot_aii, minus_ajk, aji_aik;

	// initialize
	mpfr_init2(absmax_aji, mpfr_get_prec(A[0]));
	mpfr_init2(abs_aji, mpfr_get_prec(A[0]));
	mpfr_init2(pivot_aii, mpfr_get_prec(A[0]));
	mpfr_init2(abs_pivot_aii, mpfr_get_prec(A[0]));
	mpfr_init2(aji_aik, mpfr_get_prec(A[0]));

	// initialize pivot vector
	for(i = 0; i < dim; i++)
		pivot[i] = i;

	// A decomposition
	for(i = 0; i < dim; i++)
	{
		// partial pivoting
		//absmax_aji = abs(A[ZERO_INDEX(pivot[i], i, dim)]);
		mpfr_abs(absmax_aji, A[ZERO_INDEX(pivot[i], i, dim)], _tk_default_rmode);

		max_j = i;
		for(j = i + 1; j < dim; j++)
		{
			//abs_aji = A[ZERO_INDEX(pivot[j], i, dim)];
			mpfr_abs(abs_aji, A[ZERO_INDEX(pivot[j], i, dim)], _tk_default_rmode);
			//if(absmax_aji < abs_aji)
			if(mpfr_cmp(absmax_aji, abs_aji) < 0)
			{
				max_j = j;
				//absmax_aji = abs_aji;
				mpfr_abs(absmax_aji, abs_aji, _tk_default_rmode);
			}
		}
		if(max_j != i)
		{
			tmp_index = pivot[max_j];
			pivot[max_j] = pivot[i];
			pivot[i] = tmp_index;
		}

		// select pivoted column
		//pivot_aii = A[ZERO_INDEX(pivot[i], i, dim)];
		mpfr_set(pivot_aii, A[ZERO_INDEX(pivot[i], i, dim)], _tk_default_rmode);

		// error
		//if(abs(pivot_aii) <= (double)0)
		mpfr_abs(abs_pivot_aii, pivot_aii, _tk_default_rmode);
		if(mpfr_cmp_ui(pivot_aii, 0UL))
		{
			fprintf(stderr, "mpfr_LU error!\n");

			// clear
			mpfr_clear(absmax_aji);
			mpfr_clear(abs_aji);
			mpfr_clear(pivot_aii);
			mpfr_clear(abs_pivot_aii);
			mpfr_clear(aji_aik);

			return -1;
		}

		for(j = i + 1; j < dim; j++)
		{
			//A[ZERO_INDEX(pivot[j], i, dim)] /= pivot_aii;
			mpfr_div(A[ZERO_INDEX(pivot[j], i, dim)],  A[ZERO_INDEX(pivot[j], i, dim)], pivot_aii, _tk_default_rmode);

			for(k = i + 1; k < dim; k++)
			{
				//A[ZERO_INDEX(pivot[j], k, dim)] -= A[ZERO_INDEX(pivot[j], i, dim)] * A[ZERO_INDEX(pivot[i], k, dim)];
				mpfr_mul(aji_aik, A[ZERO_INDEX(pivot[j], i, dim)], A[ZERO_INDEX(pivot[i], k, dim)], _tk_default_rmode);
				mpfr_sub(A[ZERO_INDEX(pivot[j], k, dim)],  A[ZERO_INDEX(pivot[j], k, dim)], aji_aik, _tk_default_rmode);
			}
		}
	}

	// clear
	mpfr_clear(absmax_aji);
	mpfr_clear(abs_aji);
	mpfr_clear(pivot_aii);
	mpfr_clear(abs_pivot_aii);
	mpfr_clear(aji_aik);

	return 0;

}

// rowwise only
// solve LU * x = b in x -> b := x
int mpfr_solve_LU_linear_eq(mpfr_t x[], mpfr_t LU[], mpfr_t b[], int dim, int pivot[])
{
	int i, j;
	mpfr_t luij_xj;

	mpfr_init2(luij_xj, mpfr_get_prec(x[0]));

	// x := b
	for(i = 0; i < dim; i++)
	{
		//x[i] = b[pivot[i]];
		mpfr_set(x[i], b[pivot[i]], _tk_default_rmode);
	}

	// forward substitution
	for(j = 0; j < dim; j++)
	{
		for(i = j + 1; i < dim; i++)
		{
			//x[i] -= LU[ZERO_INDEX(pivot[i], j, dim)] * x[j];
			mpfr_mul(luij_xj, LU[ZERO_INDEX(pivot[i], j, dim)], x[j], _tk_default_rmode);
			mpfr_sub(x[i], x[i], luij_xj, _tk_default_rmode);
		}
	}

	// backward substitution
	for(i = dim - 1; i >= 0; i--)
	{
		for(j = i + 1; j < dim; j++)
		{
			//x[i] -= LU[ZERO_INDEX(pivot[i], j, dim)] * x[j];
			mpfr_mul(luij_xj, LU[ZERO_INDEX(pivot[i], j, dim)], x[j], _tk_default_rmode);
			mpfr_sub(x[i], x[i], luij_xj, _tk_default_rmode);
		}

		//x[i] /= LU[ZERO_INDEX(pivot[i], i, dim)];
		mpfr_div(x[i], x[i], LU[ZERO_INDEX(pivot[i], i, dim)], _tk_default_rmode);
	}

	mpfr_clear(luij_xj);

	return 0;
}

#endif // __MPFR_H

#endif // __LINEAR_TK_C_MP__
