//******************************************************************************
// matmul_winograd.cpp : Matrix multiprication based on Winograd algorithm
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
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include "matmul_block.h"
#include "get_secv.h"

using namespace std;

// Padding to 2-powered dimensional matrix
DMatrix init_static_padding_dmatrix_winograd(DMatrix orig_mat)
{
	DMatrix ret = NULL;
	long int ret_row_dim, ret_col_dim, min_dim, i, j;

	if(orig_mat == NULL)
	{
		fprintf(stderr, "Warning: orig_mat is empty!(padding_dmatrix_winograd)\n");
		return NULL;
	}

	ret_row_dim = (long int)pow(2.0, ceil(mylog2((double)orig_mat->row_dim)));
	ret_col_dim = (long int)pow(2.0, ceil(mylog2((double)orig_mat->col_dim)));

	//printf("Padding: row_dim %ld -> %ld, col_dim %ld -> %ld\n", orig_mat->row_dim, ret_row_dim, orig_mat->col_dim, ret_col_dim);

	ret = init_dmatrix(ret_row_dim, ret_col_dim);
	if(ret == NULL)
	{
		fprintf(stderr, "Warning: padding matrix cannot be allocated!(padding_dmatrix_winograd)\n");
		return NULL;
	}

	// ret = [ A | 0 ]
	//       [---+---]
	//       [ 0 | I ]
	for(i = 0; i < orig_mat->row_dim; i++)
	{
		for(j = 0; j < orig_mat->col_dim; j++)
			set_dmatrix_ij(ret, i, j, get_dmatrix_ij(orig_mat, i, j));
	}

	min_dim = (ret_row_dim < ret_col_dim) ? ret_row_dim : ret_col_dim;
	for(i = orig_mat->row_dim; i < min_dim; i++)
		set_dmatrix_ij(ret, i, i, 1.0);

	return ret;
}

// Dynamic Padding to even dimensional matrix
DMatrix init_dynamic_padding_dmatrix_winograd(DMatrix orig_mat)
{
	DMatrix ret = NULL;
	long int ret_row_dim, ret_col_dim, min_dim, i, j;

	if(orig_mat == NULL)
	{
		fprintf(stderr, "Warning: orig_mat is empty!(padding_dmatrix_winograd)\n");
		return NULL;
	}

	if((orig_mat->row_dim % 2) >= 1)
		ret_row_dim = (long int)(orig_mat->row_dim + 1);
	if((orig_mat->col_dim % 2) >= 1)
		ret_col_dim = (long int)(orig_mat->col_dim + 1);

	//printf("Padding: row_dim %ld -> %ld, col_dim %ld -> %ld\n", orig_mat->row_dim, ret_row_dim, orig_mat->col_dim, ret_col_dim);

	ret = init_dmatrix(ret_row_dim, ret_col_dim);
	if(ret == NULL)
	{
		fprintf(stderr, "Warning: padding matrix cannot be allocated!(padding_dmatrix_winograd)\n");
		return NULL;
	}

	// ret = [ A | 0 ]
	//       [---+---]
	//       [ 0 | 1 ]
	for(i = 0; i < orig_mat->row_dim; i++)
	{
		for(j = 0; j < orig_mat->col_dim; j++)
			set_dmatrix_ij(ret, i, j, get_dmatrix_ij(orig_mat, i, j));
	}

	min_dim = (ret_row_dim < ret_col_dim) ? ret_row_dim : ret_col_dim;
	for(i = orig_mat->row_dim; i < min_dim; i++)
		set_dmatrix_ij(ret, i, i, 1.0);

	return ret;
}

// Strassen's Algorithm with static padding
void mul_dmatrix_winograd_odd_padding(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
{
	long int tmp_ret_index[4], ret_index[4];
	DMatrix tmp_ret, tmp_mat_a, tmp_mat_b;

	// padding
#ifdef STATIC_PADDING
	tmp_ret = init_static_padding_dmatrix_winograd(ret);
	tmp_mat_a = init_static_padding_dmatrix_winograd(mat_a);
	tmp_mat_b = init_static_padding_dmatrix_winograd(mat_b);
#else
	tmp_ret = init_dynamic_padding_dmatrix_winograd(ret);
	tmp_mat_a = init_dynamic_padding_dmatrix_winograd(mat_a);
	tmp_mat_b = init_dynamic_padding_dmatrix_winograd(mat_b);
#endif

	// winograd
	mul_dmatrix_winograd_even(tmp_ret, tmp_mat_a, tmp_mat_b, min_dim);

	// substitute
	tmp_ret_index[0] = 0;
	tmp_ret_index[1] = ret->row_dim;
	tmp_ret_index[2] = 0;
	tmp_ret_index[3] = ret->col_dim;
	ret_index[0] = 0;
	ret_index[1] = ret->row_dim;
	ret_index[2] = 0;
	ret_index[3] = ret->col_dim;

	subst_dmatrix_partial(ret, ret_index, tmp_ret, tmp_ret_index);

	// free
	free_dmatrix(tmp_ret);
	free_dmatrix(tmp_mat_a);
	free_dmatrix(tmp_mat_b);
}

// Fit dimension to be multiple of min_dim
void mul_dmatrix_winograd(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
{
	long int dim;

	// initialize
	dim = ret->row_dim;

	// normal matrix multiplication in case of ret_dim <= 4 
	// normal matrix multiplication in case of ret_dim <= 4 
	if((ret->row_dim <= min_dim) && (ret->col_dim <= min_dim))
	{
#ifdef USE_IMKL
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, 1.0, mat_a->element, dim, mat_b->element, dim, 0.0, ret->element, dim);
#else
		mul_dmatrix(ret, mat_a, mat_b);
#endif
		return;
	}

	// dynamic peeling in case of odd dim
	// [ A11   a12 ] [ B11   b12 ] = [ A11*B11 + a12 * b21^T   A11*b12 + a12 * b22    ]
	// [ a21^T a22 ] [ b21^T b22 ]   [ a21^T*B11 + a22 * b21^T a21^T * b12 + a22 * b22]
	if((ret->row_dim % 2 == 1) || (ret->col_dim % 2 == 1))
	{
		//printf("%d is odd -> ", ret->row_dim);
		mul_dmatrix_winograd_odd_peeling(ret, mat_a, mat_b, min_dim);
		//printf("end\n");
	}
	// normal winograd algorithm in case of even dim
	else
	{
		//printf("%d is even -> ", ret->row_dim);
		mul_dmatrix_winograd_even(ret, mat_a, mat_b, min_dim);
		//printf("end\n");
	}
}

// Strassen's Algorithm
void mul_dmatrix_winograd_odd_peeling(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
{
	long int i, j, dim, dim_h;
	DMatrix mat_a11, mat_b11, mat_c11, mat_tmp;
	DVector vec_a12, vec_a21, vec_b12, vec_b21, vec_c12, vec_c21, vec_tmp;
	double a22, b22, c22, tmp;

	// initialize
	dim = ret->row_dim;
	dim_h = dim - 1;

	if((ret->row_dim <= min_dim) && (ret->col_dim <= min_dim))
	{
		mul_dmatrix(ret, mat_a, mat_b);
		return;
	}

	// Initialize
	mat_a11 = init_dmatrix(dim_h, dim_h);
	mat_b11 = init_dmatrix(dim_h, dim_h);
	mat_c11 = init_dmatrix(dim_h, dim_h);
	mat_tmp = init_dmatrix(dim_h, dim_h);

	vec_a12 = init_dvector(dim_h);
	vec_a21 = init_dvector(dim_h);
	vec_b12 = init_dvector(dim_h);
	vec_b21 = init_dvector(dim_h);
	vec_c12 = init_dvector(dim_h);
	vec_c21 = init_dvector(dim_h);
	vec_tmp = init_dvector(dim_h);

	// set matrix elements to mat_a11, mat_b11, vec_a12, vec_b12, vec_a21, vec_b21, a22, b22
	for(i = 0; i < dim_h; i++)
	{
		set_dvector_i(vec_a12, i, get_dmatrix_ij(mat_a, i, mat_a->col_dim - 1));
		set_dvector_i(vec_b12, i, get_dmatrix_ij(mat_b, i, mat_b->col_dim - 1));
		set_dvector_i(vec_a21, i, get_dmatrix_ij(mat_a, mat_a->row_dim - 1, i));
		set_dvector_i(vec_b21, i, get_dmatrix_ij(mat_b, mat_b->row_dim - 1, i));
		for(j = 0; j < dim_h; j++)
		{
			set_dmatrix_ij(mat_a11, i, j, get_dmatrix_ij(mat_a, i, j));
			set_dmatrix_ij(mat_b11, i, j, get_dmatrix_ij(mat_b, i, j));
		}
	}
	a22 = get_dmatrix_ij(mat_a, mat_a->row_dim - 1, mat_a->col_dim - 1);
	b22 = get_dmatrix_ij(mat_b, mat_b->row_dim - 1, mat_b->col_dim - 1);

	// dynamic peeling in case of odd dim
	// [ A11   a12 ] [ B11   b12 ] = [ A11*B11 + a12 * b21^T   A11*b12 + a12 * b22     ]
	// [ a21^T a22 ] [ b21^T b22 ]   [ a21^T*B11 + a22 * b21^T a21^T * b12 + a22 * b22 ]

	// A11 * B11
	mul_dmatrix_winograd(mat_c11, mat_a11, mat_b11, min_dim);

	// a12 * b21^T
	for(i = 0; i < dim_h; i++)
	{
		for(j = 0; j < dim_h; j++)
		{
;			set_dmatrix_ij(mat_tmp, i, j, get_dvector_i(vec_a12, i) * get_dvector_i(vec_b21, j));
		}
	}
	
	add_dmatrix(mat_c11, mat_c11, mat_tmp);
	for(i = 0; i < dim_h; i++)
		for(j = 0; j < dim_h; j++)
			set_dmatrix_ij(ret, i, j, get_dmatrix_ij(mat_c11, i, j));

	// A11 * b12 + b22 * a12
	mul_dmatrix_dvec(vec_c12, mat_a11, vec_b12);
	cmul_dvector(vec_tmp, b22, vec_a12);
	add_dvector(vec_c12, vec_c12, vec_tmp);
	for(i = 0; i < dim_h; i++)
		set_dmatrix_ij(ret, i, ret->col_dim - 1, get_dvector_i(vec_c12, i));

	// a21^T*B11 + a22 * b21^T
	mul_dmatrixt_dvec(vec_c21, mat_b11, vec_a21);
	cmul_dvector(vec_tmp, a22, vec_b21);
	add_dvector(vec_c21, vec_c21, vec_tmp);
	for(i = 0; i < dim_h; i++)
		set_dmatrix_ij(ret, ret->row_dim - 1, i, get_dvector_i(vec_c21, i));

	// a21^T * b12 + a22 * b22
	c22 = ip_dvector(vec_a21, vec_b12);
	c22 += a22 * b22;
	set_dmatrix_ij(ret, ret->row_dim - 1, ret->col_dim - 1, c22);

	// free
	free_dmatrix(mat_a11);
	free_dmatrix(mat_b11);
	free_dmatrix(mat_c11);
	free_dmatrix(mat_tmp);

	free_dvector(vec_a12);
	free_dvector(vec_a21);
	free_dvector(vec_b12);
	free_dvector(vec_b21);
	free_dvector(vec_c12);
	free_dvector(vec_c21);
	free_dvector(vec_tmp);
}

// Winograd Variant of Strassen's Algorithm
void mul_dmatrix_winograd_even(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
{
//	long int min_dim = 4; // = 2^2
	DMatrix mat_s[8], mat_m[7], mat_t[2], mat_tmp_a[4], mat_tmp_b[4], mat_tmp_c[4];
	long int row_dim_h, row_dim, col_dim_h, col_dim, mid_dim, mid_dim_h;
	long int ret_index[4], mat_tmp_a_index[4], mat_tmp_b_index[4];
	long int mat_a_index[4][4], mat_b_index[4][4], mat_c_index[4][4];
	long int i;

	// initialize
	row_dim = ret->row_dim;
	col_dim = ret->col_dim;
	mid_dim = mat_a->col_dim;
	if(mid_dim != mat_b->row_dim)
	{
		fprintf(stderr, "ERROR: mat_a's col_dim %ld does not just fit mat_b's row_dim %ld!(mul_dmatrix_winograd_even)\n", mat_a->col_dim, mat_b->row_dim);
		return;
	}

	row_dim_h = row_dim / 2;
	col_dim_h = col_dim / 2;
	mid_dim_h = mid_dim / 2;

	// normal matrix multiplication in case of ret_dim <= 4 
//	if((ret->row_dim <= min_dim) && (ret->col_dim <= min_dim) && (mid_dim <= min_dim))
	if((ret->row_dim <= min_dim) || (ret->col_dim <= min_dim) || (mid_dim <= min_dim))
	{
#ifdef USE_IMKL
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row_dim, col_dim, mid_dim, 1.0, mat_a->element, mat_a->row_dim, mat_b->element, mat_b->row_dim, 0.0, ret->element, ret->row_dim);
#else
		mul_dmatrix(ret, mat_a, mat_b);
#endif
		return;
	}

	for(i = 0; i < 4; i++)
	{
		mat_s[i] = init_dmatrix(row_dim_h, mid_dim_h);
		mat_tmp_a[i] = init_dmatrix(row_dim_h, mid_dim_h);

		mat_s[i+4] = init_dmatrix(mid_dim_h, col_dim_h);
		mat_tmp_b[i] = init_dmatrix(mid_dim_h, col_dim_h);

		mat_tmp_c[i] = init_dmatrix(row_dim_h, col_dim_h);
	}
	for(i = 0; i < 7; i++)
		mat_m[i] = init_dmatrix(row_dim_h, col_dim_h);

	mat_t[0] = init_dmatrix(row_dim_h, col_dim_h);
	mat_t[1] = init_dmatrix(row_dim_h, col_dim_h);

	ret_index[0] = 0;
	ret_index[1] = row_dim_h;
	ret_index[2] = 0;
	ret_index[3] = col_dim_h;

	// indeces for A
	mat_a_index[0][0] = 0        ; mat_a_index[0][1] = row_dim_h; mat_a_index[0][2] = 0        ; mat_a_index[0][3] = mid_dim_h;
	mat_a_index[1][0] = 0        ; mat_a_index[1][1] = row_dim_h; mat_a_index[1][2] = mid_dim_h; mat_a_index[1][3] = mid_dim;
	mat_a_index[2][0] = row_dim_h; mat_a_index[2][1] = row_dim  ; mat_a_index[2][2] = 0        ; mat_a_index[2][3] = mid_dim_h;
	mat_a_index[3][0] = row_dim_h; mat_a_index[3][1] = row_dim  ; mat_a_index[3][2] = mid_dim_h; mat_a_index[3][3] = mid_dim;

	// indeces for B
	mat_b_index[0][0] = 0        ; mat_b_index[0][1] = mid_dim_h; mat_b_index[0][2] = 0        ; mat_b_index[0][3] = col_dim_h;
	mat_b_index[1][0] = 0        ; mat_b_index[1][1] = mid_dim_h; mat_b_index[1][2] = col_dim_h; mat_b_index[1][3] = col_dim;
	mat_b_index[2][0] = mid_dim_h; mat_b_index[2][1] = mid_dim  ; mat_b_index[2][2] = 0        ; mat_b_index[2][3] = col_dim_h;
	mat_b_index[3][0] = mid_dim_h; mat_b_index[3][1] = mid_dim  ; mat_b_index[3][2] = col_dim_h; mat_b_index[3][3] = col_dim;

	// indeces for C
	mat_c_index[0][0] = 0        ; mat_c_index[0][1] = row_dim_h; mat_c_index[0][2] = 0        ; mat_c_index[0][3] = col_dim_h;
	mat_c_index[1][0] = 0        ; mat_c_index[1][1] = row_dim_h; mat_c_index[1][2] = col_dim_h; mat_c_index[1][3] = col_dim;
	mat_c_index[2][0] = row_dim_h; mat_c_index[2][1] = row_dim  ; mat_c_index[2][2] = 0        ; mat_c_index[2][3] = col_dim_h;
	mat_c_index[3][0] = row_dim_h; mat_c_index[3][1] = row_dim  ; mat_c_index[3][2] = col_dim_h; mat_c_index[3][3] = col_dim;

	mat_tmp_a_index[0] = 0;
	mat_tmp_a_index[1] = row_dim_h;
	mat_tmp_a_index[2] = 0;
	mat_tmp_a_index[3] = mid_dim_h;

	mat_tmp_b_index[0] = 0;
	mat_tmp_b_index[1] = mid_dim_h;
	mat_tmp_b_index[2] = 0;
	mat_tmp_b_index[3] = col_dim_h;

	// mat_tmp_a[0], [1], [2], [3] := A11, A12, A21, A22
	// mat_tmp_b[0], [1], [2], [3] := B11, B12, B21, B22
	for(i = 0; i < 4; i++)
	{
		subst_dmatrix_partial(mat_tmp_a[i], mat_tmp_a_index, mat_a, mat_a_index[i]);
		subst_dmatrix_partial(mat_tmp_b[i], mat_tmp_b_index, mat_b, mat_b_index[i]);
	}

	// -------------------------------
	// S1 := A21 + A22
	//--------------------------------
	add_dmatrix(mat_s[0], mat_tmp_a[2], mat_tmp_a[3]);

	// -------------------------------
	// S2 := S1 - A11
	//--------------------------------
	sub_dmatrix(mat_s[1], mat_s[0], mat_tmp_a[0]);

	// -------------------------------
	// S3 := A11 - A21
	//--------------------------------
	sub_dmatrix(mat_s[2], mat_tmp_a[0], mat_tmp_a[2]);

	// -------------------------------
	// S4 := A12 - S2
	//--------------------------------
	sub_dmatrix(mat_s[3], mat_tmp_a[1], mat_s[1]);

	// -------------------------------
	// S5 := B12 - B11
	//--------------------------------
	sub_dmatrix(mat_s[4], mat_tmp_b[1], mat_tmp_b[0]);

	// -------------------------------
	// S6 := B22 - S5
	//--------------------------------
	sub_dmatrix(mat_s[5], mat_tmp_b[3], mat_s[4]);

	// -------------------------------
	// S7 := B22 - B12
	//--------------------------------
	sub_dmatrix(mat_s[6], mat_tmp_b[3], mat_tmp_b[1]);

	// -------------------------------
	// S8 := S6 - B21
	//--------------------------------
	sub_dmatrix(mat_s[7], mat_s[5], mat_tmp_b[2]);

	// -------------------------------
	// M1 := S2 * S6
	// -------------------------------
	mul_dmatrix_winograd(mat_m[0], mat_s[1], mat_s[5], min_dim);

	// -------------------------------
	// M2 := A11 * B11
	// -------------------------------
	mul_dmatrix_winograd(mat_m[1], mat_tmp_a[0], mat_tmp_b[0], min_dim);

	// -------------------------------
	// M3 := A12 * B21
	// -------------------------------
	mul_dmatrix_winograd(mat_m[2], mat_tmp_a[1], mat_tmp_b[2], min_dim);

	// -------------------------------
	// M4 := S3 * S7
	// -------------------------------
	mul_dmatrix_winograd(mat_m[3], mat_s[2], mat_s[6], min_dim);

	// -------------------------------
	// M5 := S1 * S5
	// -------------------------------
	mul_dmatrix_winograd(mat_m[4], mat_s[0], mat_s[4], min_dim);

	// -------------------------------
	// M6 := S4 * B22
	// -------------------------------
	mul_dmatrix_winograd(mat_m[5], mat_s[3], mat_tmp_b[3], min_dim);

	// -------------------------------
	// M7 := A22 * S8
	// -------------------------------
	mul_dmatrix_winograd(mat_m[6], mat_tmp_a[3], mat_s[7], min_dim);

	// -------------------------------
	// T1 := M1 + M2
	// -------------------------------
	add_dmatrix(mat_t[0], mat_m[0], mat_m[1]);

	// -------------------------------
	// T2 := T1 + M4
	// -------------------------------
	add_dmatrix(mat_t[1], mat_t[0], mat_m[3]);

	// -------------------------------
	// C11 := M2 + M3
	// -------------------------------
	add_dmatrix(mat_tmp_c[0], mat_m[1], mat_m[2]);

	// -------------------------------
	// C12 := T1 + M5 + M6
	// -------------------------------
	add_dmatrix(mat_tmp_c[1], mat_t[0], mat_m[4]);
	add_dmatrix(mat_tmp_c[1], mat_tmp_c[1], mat_m[5]);

	// -------------------------------
	// C21 := T2 - M7
	// -------------------------------
	sub_dmatrix(mat_tmp_c[2], mat_t[1], mat_m[6]);

	// -------------------------------
	// C22 := T2 + M5
	// -------------------------------
	add_dmatrix(mat_tmp_c[3], mat_t[1], mat_m[4]);

	// -------------------------------
	// RET := [C11 C12]
	//        [C21 C22]
	// -------------------------------
	for(i = 0; i < 4; i++)
		subst_dmatrix_partial(ret, mat_c_index[i], mat_tmp_c[i], ret_index);

	// free
	for(i = 0; i < 4; i++)
	{
		free_dmatrix(mat_s[i]);
		free_dmatrix(mat_tmp_a[i]);
		free_dmatrix(mat_s[i + 4]);
		free_dmatrix(mat_tmp_b[i]);
		free_dmatrix(mat_tmp_c[i]);
	}
	for(i = 0; i < 7; i++)
		free_dmatrix(mat_m[i]);
	
	free_dmatrix(mat_t[0]);
	free_dmatrix(mat_t[1]);
}

int main(int argc, char *argv[])
{
	int i, j, min_dim, max_dim, dim, iter, max_iter = 10, block_size, num_threads;
	DMatrix mat_a, mat_b, mat_c;
	double stime, etime, mat_c_normf;

	if(argc < 3)
	{
		cout << "Usage: " << argv[0] << " [min. dimension]  [max.dimension]  [block_size]"<< endl;
		return 0;
	}

	min_dim = atoi(argv[1]);
	max_dim = atoi(argv[2]);
	block_size = atoi(argv[3]);

	if(min_dim <= 0)
	{
		cout << "Illegal dimension! (min_dim = " << min_dim << ")" << endl;
		return 1;
}

	cout << "Use WINOGRAD algorithm" << endl;
	cout << setw(5) << "  dim :     SECONDS Mat.KB ||C||_F" << endl;

	// main loop
	//for(dim = min_dim; dim <= max_dim; dim += 128)
	for(dim = min_dim; dim <= max_dim; dim += 16)
	{

		// initialize matrices
		//mat_a = new double[dim * dim];
		mat_a = init_dmatrix(dim, dim);
		mat_b = init_dmatrix(dim, dim);
		mat_c = init_dmatrix(dim, dim);

		// set mat_a and mat_b
		for(i = 0; i < dim; i++)
		{
			for(j = 0; j < dim; j++)
			{
				mat_a->element[i * dim + j] = sqrt(5.0) * (double)(i + j + 1);
				mat_b->element[i * dim + j] = sqrt(3.0) * (double)(dim - (i + j));
			//	mat_a->element[i * dim + j] = 1.0 / (double)(i + j + 1);
			//	mat_b->element[i * dim + j] = (double)(i + j + 1);
			}
		}

		// matrix multiplication
		max_iter = 3;
		do
		{
			stime = get_real_secv();

			for(iter = 0; iter < max_iter; iter++)
				mul_dmatrix_winograd(mat_c, mat_a, mat_b, block_size);

			etime = get_real_secv(); etime -= stime;

			if(etime >= 1.0) break;
			max_iter *= 2;
		} while(0);

		etime /= (double)max_iter;
		mat_c_normf = normf_dmatrix(mat_c);
 
		// output
		cout << setw(5) << dim << " : " << setw(10) << setprecision(5) << etime << " " << byte_double_sqmat(dim) / 1024 << " " << mat_c_normf << endl;

		/*
		for(i = 0; i < dim; i++)
		{
			cout << " [ ";
			for(j = 0; j < dim; j++)
				cout << scientific << setprecision(3) << setw(10) << mat_a->element[i * dim + j] << " ";
			cout << " ] [ " << scientific << setprecision(3) << setw(10) << vec_x[i] << " ]   [ " << vec_b[i] << endl;
		}
		*/

		// delete matrices
		//delete mat_a;
		free_dmatrix(mat_a);
		free_dmatrix(mat_b);
		free_dmatrix(mat_c);

	} // end of main loop

	return 0;
}
