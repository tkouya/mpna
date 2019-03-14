/*************************************************/
/* matmul_block.cc : 実行列×実行列
/* [Intel] icpc matmul.cc
/* [GCC  ] g++ matmul.cc
/*************************************************/
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "matmul_block.h"
#include "get_secv.h"

using namespace std;

// Padding to 2-powered dimensional matrix
DMatrix init_static_padding_dmatrix_strassen(DMatrix orig_mat)
{
	DMatrix ret = NULL;
	long int ret_row_dim, ret_col_dim, min_dim, i, j;

	if(orig_mat == NULL)
	{
		fprintf(stderr, "Warning: orig_mat is empty!(padding_dmatrix_strassen)\n");
		return NULL;
	}

	ret_row_dim = (long int)pow(2.0, ceil(mylog2((double)orig_mat->row_dim)));
	ret_col_dim = (long int)pow(2.0, ceil(mylog2((double)orig_mat->col_dim)));

	//printf("Padding: row_dim %ld -> %ld, col_dim %ld -> %ld\n", orig_mat->row_dim, ret_row_dim, orig_mat->col_dim, ret_col_dim);

	ret = init_dmatrix(ret_row_dim, ret_col_dim);
	if(ret == NULL)
	{
		fprintf(stderr, "Warning: padding matrix cannot be allocated!(padding_dmatrix_strassen)\n");
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
DMatrix init_dynamic_padding_dmatrix_strassen(DMatrix orig_mat)
{
	DMatrix ret = NULL;
	long int ret_row_dim, ret_col_dim, min_dim, i, j;

	if(orig_mat == NULL)
	{
		fprintf(stderr, "Warning: orig_mat is empty!(padding_dmatrix_strassen)\n");
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
		fprintf(stderr, "Warning: padding matrix cannot be allocated!(padding_dmatrix_strassen)\n");
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
void mul_dmatrix_strassen_odd_padding(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
{
	long int tmp_ret_index[4], ret_index[4];
	DMatrix tmp_ret, tmp_mat_a, tmp_mat_b;

	// padding
#ifdef STATIC_PADDING
	tmp_ret = init_static_padding_dmatrix_strassen(ret);
	tmp_mat_a = init_static_padding_dmatrix_strassen(mat_a);
	tmp_mat_b = init_static_padding_dmatrix_strassen(mat_b);
#else
	tmp_ret = init_dynamic_padding_dmatrix_strassen(ret);
	tmp_mat_a = init_dynamic_padding_dmatrix_strassen(mat_a);
	tmp_mat_b = init_dynamic_padding_dmatrix_strassen(mat_b);
#endif

	// strassen
	mul_dmatrix_strassen_even(tmp_ret, tmp_mat_a, tmp_mat_b, min_dim);

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
void mul_dmatrix_strassen(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
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
		mul_dmatrix_strassen_odd_peeling(ret, mat_a, mat_b, min_dim);
		//printf("end\n");
	}
	// normal strassen algorithm in case of even dim
	else
	{
		//printf("%d is even -> ", ret->row_dim);
		mul_dmatrix_strassen_even(ret, mat_a, mat_b, min_dim);
		//printf("end\n");
	}
}

// Strassen's Algorithm
void mul_dmatrix_strassen_odd_peeling(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
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
	mul_dmatrix_strassen(mat_c11, mat_a11, mat_b11, min_dim);

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


// Strassen's Algorithm
void mul_dmatrix_strassen_even(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
{
//	long int min_dim = 4; // = 2^2
	DMatrix mat_p[7], mat_tmp_a[7], mat_tmp_b[7], mat_tmp_c[4];
	long int dim_h, dim;
	long int ret_index[4], mat_ar_index[7][4], mat_al_index[7][4], mat_br_index[7][4], mat_bl_index[7][4], mat_c_index[4][4];
	long int i;

	// initialize
	dim = ret->row_dim;
	dim_h = dim / 2;

	// normal matrix multiplication in case of ret_dim <= 4 
	if((ret->row_dim <= min_dim) && (ret->col_dim <= min_dim))
	{
		mul_dmatrix(ret, mat_a, mat_b);
		return;
	}

	for(i = 0; i < 7; i++)
	{
		mat_p[i] = init_dmatrix(dim_h, dim_h);
		mat_tmp_a[i] = init_dmatrix(dim_h, dim_h);
		mat_tmp_b[i] = init_dmatrix(dim_h, dim_h);
	}
	for(i = 0; i < 4; i++)
		mat_tmp_c[i] = init_dmatrix(dim_h, dim_h);

	ret_index[0] = 0;
	ret_index[1] = dim_h;
	ret_index[2] = 0;
	ret_index[3] = dim_h;

	// -------------------------------
	// P1 := (A11 + A22) * (B11 + B22)
	//--------------------------------

	// A11 + A22
	mat_ar_index[0][0] = 0;
	mat_ar_index[0][1] = dim_h;
	mat_ar_index[0][2] = 0;
	mat_ar_index[0][3] = dim_h;

	mat_al_index[0][0] = dim_h;
	mat_al_index[0][1] = dim;
	mat_al_index[0][2] = dim_h;
	mat_al_index[0][3] = dim;
	//printf("P1a:\n");
	add_dmatrix_partial(mat_tmp_a[0], ret_index, mat_a, mat_ar_index[0], mat_a, mat_al_index[0]);

	// B11 + B22
	mat_br_index[0][0] = 0;
	mat_br_index[0][1] = dim_h;
	mat_br_index[0][2] = 0;
	mat_br_index[0][3] = dim_h;

	mat_bl_index[0][0] = dim_h;
	mat_bl_index[0][1] = dim;
	mat_bl_index[0][2] = dim_h;
	mat_bl_index[0][3] = dim;
	//printf("P1b:\n");
	add_dmatrix_partial(mat_tmp_b[0], ret_index, mat_b, mat_br_index[0], mat_b, mat_bl_index[0]);

	// P1 = tmp_a * tmp_b
	//printf("P1:\n");
	mul_dmatrix_strassen(mat_p[0], mat_tmp_a[0], mat_tmp_b[0], min_dim);

	// -------------------------------
	// P2 := (A21 + A22) * B11
	// -------------------------------

	// A21 + A22
	mat_ar_index[1][0] = dim_h;
	mat_ar_index[1][1] = dim;
	mat_ar_index[1][2] = 0;
	mat_ar_index[1][3] = dim_h;

	mat_al_index[1][0] = dim_h;
	mat_al_index[1][1] = dim;
	mat_al_index[1][2] = dim_h;
	mat_al_index[1][3] = dim;
	add_dmatrix_partial(mat_tmp_a[1], ret_index, mat_a, mat_ar_index[1], mat_a, mat_al_index[1]);

	// B11
	mat_br_index[1][0] = 0;
	mat_br_index[1][1] = dim_h;
	mat_br_index[1][2] = 0;
	mat_br_index[1][3] = dim_h;
	subst_dmatrix_partial(mat_tmp_b[1], ret_index, mat_b, mat_br_index[1]);

	// P2 = tmp_a * tmp_b
	//printf("P2:\n");
	mul_dmatrix_strassen(mat_p[1], mat_tmp_a[1], mat_tmp_b[1], min_dim);

	// -------------------------------
	// P3 := A11 * (B12 - B22)
	// -------------------------------

	// A11
	mat_ar_index[2][0] = 0;
	mat_ar_index[2][1] = dim_h;
	mat_ar_index[2][2] = 0;
	mat_ar_index[2][3] = dim_h;
	subst_dmatrix_partial(mat_tmp_a[2], ret_index, mat_a, mat_ar_index[2]);

	// B12 - B22
	mat_br_index[2][0] = 0;
	mat_br_index[2][1] = dim_h;
	mat_br_index[2][2] = dim_h;
	mat_br_index[2][3] = dim;

	mat_bl_index[2][0] = dim_h;
	mat_bl_index[2][1] = dim;
	mat_bl_index[2][2] = dim_h;
	mat_bl_index[2][3] = dim;
	sub_dmatrix_partial(mat_tmp_b[2], ret_index, mat_b, mat_br_index[2], mat_b, mat_bl_index[2]);

	// P3 = tmp_a * tmp_b
	//printf("P3:\n");
	mul_dmatrix_strassen(mat_p[2], mat_tmp_a[2], mat_tmp_b[2], min_dim);

	// -------------------------------
	// P4 := A22 * (B21 - B11)
	// -------------------------------

	// A22
	mat_ar_index[3][0] = dim_h;
	mat_ar_index[3][1] = dim;
	mat_ar_index[3][2] = dim_h;
	mat_ar_index[3][3] = dim;
	subst_dmatrix_partial(mat_tmp_a[3], ret_index, mat_a, mat_ar_index[3]);

	// B21 - B11
	mat_br_index[3][0] = dim_h;
	mat_br_index[3][1] = dim;
	mat_br_index[3][2] = 0;
	mat_br_index[3][3] = dim_h;

	mat_bl_index[3][0] = 0;
	mat_bl_index[3][1] = dim_h;
	mat_bl_index[3][2] = 0;
	mat_bl_index[3][3] = dim_h;
	sub_dmatrix_partial(mat_tmp_b[3], ret_index, mat_b, mat_br_index[3], mat_b, mat_bl_index[3]);

	// P4 = tmp_a * tmp_b
	//printf("P4:\n");
	mul_dmatrix_strassen(mat_p[3], mat_tmp_a[3], mat_tmp_b[3], min_dim);

	// -------------------------------
	// P5 := (A11 + A12) * B22
	// -------------------------------

	// A11 + A12
	mat_ar_index[4][0] = 0;
	mat_ar_index[4][1] = dim_h;
	mat_ar_index[4][2] = 0;
	mat_ar_index[4][3] = dim_h;

	mat_al_index[4][0] = 0;
	mat_al_index[4][1] = dim_h;
	mat_al_index[4][2] = dim_h;
	mat_al_index[4][3] = dim;
	add_dmatrix_partial(mat_tmp_a[4], ret_index, mat_a, mat_ar_index[4], mat_a, mat_al_index[4]);

	// B22
	mat_br_index[4][0] = dim_h;
	mat_br_index[4][1] = dim;
	mat_br_index[4][2] = dim_h;
	mat_br_index[4][3] = dim;
	subst_dmatrix_partial(mat_tmp_b[4], ret_index, mat_b, mat_br_index[4]);

	// P5 = tmp_a * tmp_b
	//printf("P5:\n");
	mul_dmatrix_strassen(mat_p[4], mat_tmp_a[4], mat_tmp_b[4], min_dim);

	// -------------------------------
	// P6 := (A21 - A11) * (B11 + B12)
	// -------------------------------
	// A21 - A11
	mat_ar_index[5][0] = dim_h;
	mat_ar_index[5][1] = dim;
	mat_ar_index[5][2] = 0;
	mat_ar_index[5][3] = dim_h;

	mat_al_index[5][0] = 0;
	mat_al_index[5][1] = dim_h;
	mat_al_index[5][2] = 0;
	mat_al_index[5][3] = dim_h;
	sub_dmatrix_partial(mat_tmp_a[5], ret_index, mat_a, mat_ar_index[5], mat_a, mat_al_index[5]);

	// B11 + B12
	mat_br_index[5][0] = 0;
	mat_br_index[5][1] = dim_h;
	mat_br_index[5][2] = 0;
	mat_br_index[5][3] = dim_h;

	mat_bl_index[5][0] = 0;
	mat_bl_index[5][1] = dim_h;
	mat_bl_index[5][2] = dim_h;
	mat_bl_index[5][3] = dim;
	add_dmatrix_partial(mat_tmp_b[5], ret_index, mat_b, mat_br_index[5], mat_b, mat_bl_index[5]);

	// P6 = tmp_a * tmp_b
	//printf("P6:\n");
	mul_dmatrix_strassen(mat_p[5], mat_tmp_a[5], mat_tmp_b[5], min_dim);

	// -------------------------------
	// P7 := (A12 - A22) * (B21 + B22)
	// -------------------------------

	// A12 - A22
	mat_ar_index[6][0] = 0;
	mat_ar_index[6][1] = dim_h;
	mat_ar_index[6][2] = dim_h;
	mat_ar_index[6][3] = dim;

	mat_al_index[6][0] = dim_h;
	mat_al_index[6][1] = dim;
	mat_al_index[6][2] = dim_h;
	mat_al_index[6][3] = dim;
	sub_dmatrix_partial(mat_tmp_a[6], ret_index, mat_a, mat_ar_index[6], mat_a, mat_al_index[6]);

	// B21 + B22
	mat_br_index[6][0] = dim_h;
	mat_br_index[6][1] = dim;
	mat_br_index[6][2] = 0;
	mat_br_index[6][3] = dim_h;

	mat_bl_index[6][0] = dim_h;
	mat_bl_index[6][1] = dim;
	mat_bl_index[6][2] = dim_h;
	mat_bl_index[6][3] = dim;
	add_dmatrix_partial(mat_tmp_b[6], ret_index, mat_b, mat_br_index[6], mat_b, mat_bl_index[6]);

	// P7 = tmp_a * tmp_b
	//printf("P7:\n");
	mul_dmatrix_strassen(mat_p[6], mat_tmp_a[6], mat_tmp_b[6], min_dim);

	// -------------------------------
	// C11 := P1 + P4 - P5 + P7
	// -------------------------------
	add_dmatrix(mat_tmp_c[0], mat_p[0], mat_p[3]);
	sub_dmatrix(mat_tmp_c[0], mat_tmp_c[0], mat_p[4]);
	add_dmatrix(mat_tmp_c[0], mat_tmp_c[0], mat_p[6]);
	mat_c_index[0][0] = 0;
	mat_c_index[0][1] = dim_h;
	mat_c_index[0][2] = 0;
	mat_c_index[0][3] = dim_h;
	subst_dmatrix_partial(ret, mat_c_index[0], mat_tmp_c[0], ret_index);

	// -------------------------------
	// C12 := P3 + P5
	// -------------------------------
	add_dmatrix(mat_tmp_c[1], mat_p[2], mat_p[4]);
	mat_c_index[1][0] = 0;
	mat_c_index[1][1] = dim_h;
	mat_c_index[1][2] = dim_h;
	mat_c_index[1][3] = dim;
	subst_dmatrix_partial(ret, mat_c_index[1], mat_tmp_c[1], ret_index);

	// -------------------------------
	// C21 := P2 + P4
	// -------------------------------
	add_dmatrix(mat_tmp_c[2], mat_p[1], mat_p[3]);
	mat_c_index[2][0] = dim_h;
	mat_c_index[2][1] = dim;
	mat_c_index[2][2] = 0;
	mat_c_index[2][3] = dim_h;
	subst_dmatrix_partial(ret, mat_c_index[2], mat_tmp_c[2], ret_index);

	// -------------------------------
	// C22 := P1 + P3 - P2 + P6
	// -------------------------------
	add_dmatrix(mat_tmp_c[3], mat_p[0], mat_p[2]);
	sub_dmatrix(mat_tmp_c[3], mat_tmp_c[3], mat_p[1]);
	add_dmatrix(mat_tmp_c[3], mat_tmp_c[3], mat_p[5]);
	mat_c_index[3][0] = dim_h;
	mat_c_index[3][1] = dim;
	mat_c_index[3][2] = dim_h;
	mat_c_index[3][3] = dim;
	subst_dmatrix_partial(ret, mat_c_index[3], mat_tmp_c[3], ret_index);

	// free
	for(i = 0; i < 7; i++)
	{
		free_dmatrix(mat_p[i]);
		free_dmatrix(mat_tmp_a[i]);
		free_dmatrix(mat_tmp_b[i]);
	}
	for(i = 0; i < 4; i++)
		free_dmatrix(mat_tmp_c[i]);
}

// Block matrix multiplicaiton
void mul_dmatrix_block(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim)
{
	int row_padding_flag = 0, col_padding_flag = 0, mid_padding_flag = 0;
	long i, j, k, row_dim, col_dim, mid_dim;
	long int num_div_row, num_div_col, num_div_mid, max_num_div;
	long int **mat_a_index, small_mat_a_index[4];
	long int **mat_b_index, small_mat_b_index[4];
	long int **ret_index, small_ret_index[4];
	//DMatrix small_ret[1024], small_mat_a[1024], small_mat_b[1024], small_tmp_mat;
	DMatrix *small_ret, *small_mat_a, *small_mat_b, *small_tmp_mat;

	// initialize
	row_dim = ret->row_dim;
	col_dim = ret->col_dim;
	mid_dim = mat_a->col_dim;
	if(mid_dim != mat_b->row_dim)
	{
		fprintf(stderr, "ERROR: mat_a's col_dim %ld does not just fit mat_b's row_dim %ld!(_bncomp_mul_dmatrix_block)\n", mat_a->col_dim, mat_b->row_dim);
		return;
	}

	// normal matrix multiplication in case of ret_dim <= 4 
	if((ret->row_dim <= min_dim) && (ret->col_dim <= min_dim) && (mid_dim <= min_dim))
	{
		mul_dmatrix(ret, mat_a, mat_b);
		return;
	}

	// Number of division of matrix
	num_div_row = (ret->row_dim) / min_dim;
	if((ret->row_dim % min_dim) >= 1)
	{
		row_padding_flag = 1;
		num_div_row++;
	}

	num_div_mid = mid_dim / min_dim;
	if((mid_dim % min_dim) >= 1)
	{
		mid_padding_flag = 1;
		num_div_mid++;
	}

	num_div_col = (ret->col_dim) / min_dim;
	if((ret->col_dim % min_dim) >= 1)
	{
		col_padding_flag = 1;
		num_div_col++;
	}

	max_num_div = (num_div_row > num_div_mid) ? num_div_row : num_div_mid;
	max_num_div = (max_num_div > num_div_col) ? max_num_div : num_div_col;

	mat_a_index = (long int **)calloc(max_num_div, sizeof(long int *));
	mat_b_index = (long int **)calloc(max_num_div, sizeof(long int *));
	ret_index = (long int **)calloc(max_num_div, sizeof(long int *));

	#pragma omp parallel for
	for(i = 0; i < max_num_div; i++)
	{
		mat_a_index[i] = (long int *)calloc(4, sizeof(long int));
		mat_b_index[i] = (long int *)calloc(4, sizeof(long int));
		ret_index[i] = (long int *)calloc(4, sizeof(long int));
	}

	small_ret = (DMatrix *)calloc(sizeof(DMatrix), num_div_col);
	small_mat_a = (DMatrix *)calloc(sizeof(DMatrix), num_div_mid);
	small_mat_b = (DMatrix *)calloc(sizeof(DMatrix), num_div_mid);
	small_tmp_mat = (DMatrix *)calloc(sizeof(DMatrix), num_div_mid);

	#pragma omp parallel for
	for(i = 0; i < num_div_col; i++)
		small_ret[i] = init_dmatrix(min_dim, min_dim);

	#pragma omp parallel for
	for(i = 0; i < num_div_mid; i++)
	{
		small_mat_a[i] = init_dmatrix(min_dim, min_dim);
		small_mat_b[i] = init_dmatrix(min_dim, min_dim);
		small_tmp_mat[i] = init_dmatrix(min_dim, min_dim);
	}

	// mail loop
	small_mat_a_index[0] = 0;
	small_mat_a_index[1] = min_dim;
	small_mat_a_index[2] = 0;
	small_mat_a_index[3] = min_dim;

	small_mat_b_index[0] = 0;
	small_mat_b_index[1] = min_dim;
	small_mat_b_index[2] = 0;
	small_mat_b_index[3] = min_dim;

	small_ret_index[0] = 0;
	small_ret_index[1] = min_dim;
	small_ret_index[2] = 0;
	small_ret_index[3] = min_dim;

	// mail loop
	for(i = 0; i < num_div_row; i++)
	{
		#pragma omp parallel for
		for(j = 0; j < num_div_mid; j++)
		{
			// copy matrices
			mat_a_index[j][0] = i * min_dim;
			mat_a_index[j][1] = (i + 1) * min_dim;
			mat_a_index[j][2] = j * min_dim;
			mat_a_index[j][3] = (j + 1) * min_dim;
			subst_dmatrix_partial_checked(small_mat_a[j], small_mat_a_index, mat_a, mat_a_index[j]);
		}

		for(j = 0; j < num_div_col; j++)
		{
			set0_dmatrix(small_ret[j]);

			#pragma omp parallel for
			for(k = 0; k < num_div_mid; k++)
			{
				// copy matrices
				mat_b_index[k][0] = k * min_dim;
				mat_b_index[k][1] = (k + 1) * min_dim;
				mat_b_index[k][2] = j * min_dim;
				mat_b_index[k][3] = (j + 1) * min_dim;
				subst_dmatrix_partial_checked(small_mat_b[k], small_mat_b_index, mat_b, mat_b_index[k]);
				//_bncomp_subst_dmatrix_partial_checked(small_mat_b[k], small_mat_b_index, mat_b, mat_b_index[k]);
				// ret[j] += small_mat_a[i][k] * small_mat_b[k][j];
				mul_dmatrix(small_tmp_mat[k], small_mat_a[k], small_mat_b[k]);
			}

			for(k = 0; k < num_div_mid; k++)
				add_dmatrix(small_ret[j], small_ret[j], small_tmp_mat[k]);

			ret_index[j][0] = i * min_dim;
			ret_index[j][1] = (i + 1) * min_dim;
			ret_index[j][2] = j * min_dim;
			ret_index[j][3] = (j + 1) * min_dim;
			//subst_dmatrix_partial_checked(ret, ret_index[j], small_ret[j], small_ret_index);
			subst_dmatrix_partial_checked(ret, ret_index[j], small_ret[j], small_ret_index);
		}
	}

	#pragma omp parallel for
	for(i = 0; i < max_num_div; i++)
	{
		free(mat_a_index[i]);
		free(mat_b_index[i]);
		free(ret_index[i]);
	}
	free(mat_a_index);
	free(mat_b_index);
	free(ret_index);

	#pragma omp parallel for
	for(i = 0; i < num_div_col; i++)
		free_dmatrix(small_ret[i]);

	#pragma omp parallel for
	for(i = 0; i < num_div_mid; i++)
	{
		free_dmatrix(small_mat_a[i]);
		free_dmatrix(small_mat_b[i]);
		free_dmatrix(small_tmp_mat[i]);
	}
	free(small_ret);
	free(small_mat_a);
	free(small_mat_b);
	free(small_tmp_mat);
}


int main(int argc, char *argv[])
{
	int i, j, min_dim, max_dim, dim, iter, max_iter = 10, block_size, num_threads;
	DMatrix mat_a, mat_b, mat_c;
	double stime, etime, mat_c_normf;

	// 次元数入力
//	cout << "DIM = ";
//	cin >> dim;

#if defined(BLOCK) || defined(STRASSEN)
	if(argc < 3)
	{
		cout << "Usage: " << argv[0] << " [min. dimension]  [max.dimension]  [block_size]"<< endl;
		return EXIT_SUCCESS;
	}

	min_dim = atoi(argv[1]);
	max_dim = atoi(argv[2]);
	block_size = atoi(argv[3]);
#else // (BLOCK || STRASSEN)
	if(argc < 3)
	{
		cout << "Usage: " << argv[0] << " [min. dimension]  [max.dimension] "<< endl;
		return EXIT_SUCCESS;
	}

	min_dim = atoi(argv[1]);
	max_dim = atoi(argv[2]);
#endif //(BLOCK || STRASSEN)

#ifdef _OPENMP
	cout << "num_threads: ";
	cin >> num_threads;

	omp_set_num_threads(num_threads);
#endif // _OPENMP

	if(min_dim <= 0)
	{
		cout << "Illegal dimension! (min_dim = " << min_dim << ")" << endl;
		return EXIT_FAILURE;
	}


#ifdef STRASSEN
	cout << "Use STRASSEN algorithm" << endl;
	cout << setw(5) << "  dim :     SECONDS Mat.KB" << endl;
#elif BLOCK // STRASSEN
	cout << "Use Block algorithm" << endl;
	cout << setw(5) << "  dim :     SECONDS GFLOPS Mat.KB" << endl;
#else // STRASSEN
	cout << "Use Simple algorithm" << endl;
	cout << setw(5) << "  dim :     SECONDS GFLOPS Mat.KB" << endl;
#endif // STRASSEN

	// mainloop
	//for(dim = min_dim; dim <= max_dim; dim += 128)
	for(dim = min_dim; dim <= max_dim; dim += 16)
	{

		// 変数初期化
		//mat_a = new double[dim * dim];
		mat_a = init_dmatrix(dim, dim);
		mat_b = init_dmatrix(dim, dim);
		mat_c = init_dmatrix(dim, dim);

		// mat_aとmat_bに値入力
		for(i = 0; i < dim; i++)
		{
			for(j = 0; j < dim; j++)
			{
			//	mat_a->element[i * dim + j] = sqrt(5.0) * (double)(i + j + 1);
			//	mat_b->element[i * dim + j] = sqrt(3.0) * (double)(dim - (i + j));
				mat_a->element[i * dim + j] = 1.0 / (double)(i + j + 1);
				mat_b->element[i * dim + j] = (double)(i + j + 1);
			}
		}

		// 行列×行列
		max_iter = 3;
		do
		{
			stime = get_real_secv();

#ifdef BLOCK
			for(iter = 0; iter < max_iter; iter++)
				mul_dmatrix_block(mat_c, mat_a, mat_b, block_size);
#elif STRASSEN
			for(iter = 0; iter < max_iter; iter++)
				mul_dmatrix_strassen(mat_c, mat_a, mat_b, block_size);
#else // SIMPLE
			for(iter = 0; iter < max_iter; iter++)
				mul_dmatrix(mat_c, mat_a, mat_b);
#endif // BLOCK & STRASSEN
			etime = get_real_secv(); etime -= stime;

			if(etime >= 1.0) break;
			max_iter *= 2;
		} while(0);

		etime /= (double)max_iter;
		mat_c_normf = normf_dmatrix(mat_c);
 
		// 出力
		//cout << "Dimension     : " << dim << " * " << dim << endl;
		//cout << "Comp.Time(sec): " << setprecision(3) << etime << endl;
		//cout << "Gflops        : " << setprecision(3) << matvec_mul_gflops(etime, dim) << endl;
//		cout << setw(5) << dim << " : " << setw(10) << setprecision(5) << matmul_gflops(etime, dim) << " " << byte_double_sqmat(dim) / 1024 << endl;
#ifdef STRASSEN
		cout << setw(5) << dim << " : " << setw(10) << setprecision(5) << etime << " " << byte_double_sqmat(dim) / 1024 << " " << mat_c_normf << endl;
#else // STRASSEN
		cout << setw(5) << dim << " : " << setw(10) << setprecision(5) << etime << " " << matmul_gflops(etime, dim) << " " << byte_double_sqmat(dim) / 1024 << " " << mat_c_normf << endl;
#endif // STRASSEN

		/*
		for(i = 0; i < dim; i++)
		{
			cout << " [ ";
			for(j = 0; j < dim; j++)
				cout << scientific << setprecision(3) << setw(10) << mat_a->element[i * dim + j] << " ";
			cout << " ] [ " << scientific << setprecision(3) << setw(10) << vec_x[i] << " ]   [ " << vec_b[i] << endl;
		}
		*/

		// 変数消去
		//delete mat_a;
		free_dmatrix(mat_a);
		free_dmatrix(mat_b);
		free_dmatrix(mat_c);

	} // end of mainloop

	return EXIT_SUCCESS;
}
