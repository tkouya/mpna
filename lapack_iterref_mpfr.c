/*************************************************/
/* linear_eq_dgetrf.c : 連立一次方程式求解
/* icc linear_eq.c -L/usr/local/lib -llapacke -llapack -lcblas -lrefblas -L/opt/intel/lib/intel64 -lifcore
/*************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lapacke/lapacke.h"
#include "cblas.h"

// array = {constant, constant, ..., constant}
void drepeat_constant_array(double array[], double constant, int dim)
{
	while(--dim >= 0)
		array[dim] = constant;
}

// Iterative ref.
// ret[dim] := A[dim * dim]^(-1) * b[dim]
int d_iterative_refinement(double ret[], double mat_a[], double vec_b[], int dim, int max_itimes, double rtol, double atol)
{
	lapack_int info, itimes, i;
	lapack_int *pivot;
	double *vec_z, *vec_r, *vec_x, *mat_lu;
	double norm_z, old_norm_z, norm_x, alpha = 0.1;

	// ピボット・ベクトル初期化
	pivot = (lapack_int *)calloc(dim, sizeof(lapack_int));
	mat_lu = (double *)calloc(dim * dim, sizeof(double));
	vec_z = (double *)calloc(dim, sizeof(double));
	vec_r = (double *)calloc(dim, sizeof(double));
	vec_x = (double *)calloc(dim, sizeof(double));

	// mat_lu := mat_a
	cblas_dcopy(dim * dim, mat_a, 1, mat_lu, 1);

	// x0 := 0
	// r0 := b
	drepeat_constant_array(vec_x, 0.0, dim);

	norm_x = cblas_dnrm2(dim, vec_x, 1);
	cblas_dcopy(dim, vec_b, 1, vec_r, 1);
	cblas_dcopy(dim, vec_r, 1, vec_z, 1);

	old_norm_z = cblas_dnrm2(dim, vec_z, 1);

	// LU分解
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_lu, dim, pivot);
	//printf("DGETRF info = %d\n", info);

	// ループ
	for(itimes = 0; itimes < max_itimes; itimes++)
	{
		// mat_a * z = rを解く
		// 前進，後退代入
		info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', dim, 1, mat_lu, dim, pivot, vec_z, 1);
		//printf("DGETRS info = %d\n", info);

		// x += z
		cblas_daxpy(dim, 1.0, vec_z, 1, vec_x, 1);

		// 収束判定(1)
		// ||(x + z) - x|| = ||z|| <= rtol * ||x|| + atol ?
		norm_z = cblas_dnrm2(dim, vec_z, 1);
		if(norm_z <= rtol * norm_x + atol)
			break;

		// 収束判定(2)
		// ||z|| >= alpha * ||old_z||
		if(norm_z >= alpha * old_norm_z)
			break;
	
		old_norm_z = norm_z;

		// ||x||
		norm_x = cblas_dnrm2(dim, vec_x, 1);

		// r := vec_b - mat_a * x
		cblas_dcopy(dim, vec_b, 1, vec_r, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, -1.0, mat_a, dim, vec_x, 1, 1.0, vec_r, 1);

		// z := r
		cblas_dcopy(dim, vec_r, 1, vec_z, 1);
	}

	// ret := x
	cblas_dcopy(dim, vec_x, 1, ret, 1);

	// ピボット・ベクトル消去
	free(pivot);
	free(mat_lu);
	free(vec_z);
	free(vec_r);
	free(vec_x);

	return itimes;
}

// ETF
#include "etf.c"


// Dot2
double ddot2(double vec_x[], double vec_y[], int dim)
{
	int i;
	double p, s, q, r, h;

	p = two_prod(vec_x[0], vec_y[0], &s);

	for(i = 1; i < dim; i++)
	{
		h = two_prod(vec_x[i], vec_y[i], &r);
		p = two_sum(p, h, &q);
		
		s += q + r;
	}

	return p + s;
}

// log2(x) := log10(x) / log10(2)
#define DLOG2(x) (log10((x)) / 0.30102999566398119521373889472449)

// SplitMat_A
void dsplit_mat(double high_mat[], double low_mat[], double mat[], int row_dim, int col_dim)
{
	int i, j;
	int num_digits = 53; // IEEE double prec.
	double *s;
	double mu, abs_aij, t_exp;

	// 初期化
	s = (double *)calloc(row_dim * col_dim, sizeof(double));

	// mu[i] = max_j |mat[i, j]|
	for(i = 0; i < row_dim; i++)
	{
		mu = fabs(mat[i * col_dim + 0]);
		for(j = 1; j < col_dim; j++)
		{
			abs_aij = fabs(mat[i * col_dim + j]);
			if(abs_aij > mu)
				mu = abs_aij;
		}

		// t_exp = ceil(log2(mu)) + ceil(s + log2(col_dim + 1) / 2)
		t_exp = ceil(DLOG2(mu)) + ceil(( (double)num_digits + DLOG2((double)(col_dim + 1))/2.0 ));

		// s[i, j] = 2^t_exp
		for(j = 0; j < col_dim; j++)
			s[i * col_dim + j] = pow(2.0, t_exp);
	}

	// tmp_mat := mat + s
	cblas_dcopy(row_dim * col_dim, mat, 1, high_mat, 1);
	cblas_daxpy(row_dim * col_dim, 1.0, s, 1, high_mat, 1);

	// high_mat := tmp_mat - s
	cblas_daxpy(row_dim * col_dim, -1.0, s, 1, high_mat, 1);

	// low_mat := mat - high_mat
	cblas_dcopy(row_dim * col_dim, mat, 1, low_mat, 1);
	cblas_daxpy(row_dim * col_dim, -1.0, high_mat, 1, low_mat, 1);

	// 解放
	free(s);
}

// SplitMat_B
void dsplit_mat_t(double high_mat[], double low_mat[], double mat[], int row_dim, int col_dim)
{
	int i, j;
	int num_digits = 53; // IEEE double prec.
	double *s;
	double mu, abs_aij, t_exp;

	// 初期化
	s = (double *)calloc(row_dim * col_dim, sizeof(double));

	// mu[j] = max_j |mat[i, j]|
	for(j = 0; j < col_dim; j++)
	{
		mu = fabs(mat[0 * col_dim + j]);
		for(i = 1; i < row_dim; i++)
		{
			abs_aij = fabs(mat[i * col_dim + j]);
			if(abs_aij > mu)
				mu = abs_aij;
		}

		// t_exp = ceil(log2(mu)) + ceil(s + log2(col_dim + 1) / 2)
		t_exp = ceil(DLOG2(mu)) + ceil(((double)num_digits + DLOG2((double)(row_dim + 1)) / 2.0));

		// s[i, j] = 2^t_exp
		for(i = 0; i < row_dim; i++)
			s[i * col_dim + j] = pow(2.0, t_exp);
	}

	// tmp_mat := mat + s
	cblas_dcopy(row_dim * col_dim, mat, 1, high_mat, 1);
	cblas_daxpy(row_dim * col_dim, 1.0, s, 1, high_mat, 1);

	// high_mat := tmp_mat - s
	cblas_daxpy(row_dim * col_dim, -1.0, s, 1, high_mat, 1);

	// low_mat := mat - high_mat
	cblas_dcopy(row_dim * col_dim, mat, 1, low_mat, 1);
	cblas_daxpy(row_dim * col_dim, -1.0, high_mat, 1, low_mat, 1);

	// 出力
/*	printf("high_mat: \n");
	for(i = 0; i < row_dim; i++)
	{
		for(j = 0; j < col_dim; j++)
			printf("%25.17e ", high_mat[i * col_dim + j]);
		printf("\n");
	}
	printf("low_mat: \n");
	for(i = 0; i < row_dim; i++)
	{
		for(j = 0; j < col_dim; j++)
			printf("%25.17e ", low_mat[i * col_dim + j]);
		printf("\n");
	}
	printf("high_mat + low_mat: \n");
	for(i = 0; i < row_dim; i++)
	{
		for(j = 0; j < col_dim; j++)
			printf("%25.17e ", high_mat[i * col_dim + j] + low_mat[i * col_dim + j]);
		printf("\n");
	}
*/
	// 解放
	free(s);
}


// Accurate_MM
// ret[row_dim * col_dim] := mat_a[row_dim * mid_dim] * mat_b[mid_dim * col_dim]
void accurate_dgemm(double ret[], double mat_a[], double mat_b[], int row_dim, int col_dim, int mid_dim)
{
	double *high_mat_a, *low_mat_a, *high_mat_b, *low_mat_b, *tmp_mat;

	// 初期化
	high_mat_a = (double *)calloc(row_dim * mid_dim, sizeof(double));
	low_mat_a = (double *)calloc(row_dim * mid_dim, sizeof(double));
	high_mat_b = (double *)calloc(mid_dim * col_dim, sizeof(double));
	low_mat_b = (double *)calloc(mid_dim * col_dim, sizeof(double));
	//tmp_mat = (double *)calloc(row_dim * col_dim, sizeof(double));

	// ret, tmp_mat := 0
	drepeat_constant_array(ret, 0.0, row_dim * col_dim);
	//drepeat_constant_array(tmp_mat, 0.0, row_dim * col_dim);

	// 分割
	dsplit_mat(high_mat_a, low_mat_a, mat_a, row_dim, mid_dim);
	dsplit_mat_t(high_mat_b, low_mat_b, mat_b, mid_dim, col_dim);

	// ret += high_mat_a * low_mat_b + low_mat_a * mat_b
	// ret += high_mat_a * high_mat_b
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row_dim, col_dim, mid_dim, 1.0, high_mat_a, row_dim, low_mat_b, mid_dim, 1.0, ret, row_dim);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row_dim, col_dim, mid_dim, 1.0, low_mat_a, row_dim, mat_b, mid_dim, 1.0, ret, row_dim);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row_dim, col_dim, mid_dim, 1.0, high_mat_a, row_dim, high_mat_b, mid_dim, 1.0, ret, row_dim);

	// tmp_mat := high_mat_a * low_mat_b + low_mat_a * mat_b
	// ret += tmp_mat
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row_dim, col_dim, mid_dim, 1.0, high_mat_a, row_dim, low_mat_b, mid_dim, 1.0, tmp_mat, row_dim);
	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row_dim, col_dim, mid_dim, 1.0, low_mat_a, row_dim, mat_b, mid_dim, 1.0, tmp_mat, row_dim);
	//cblas_daxpy(row_dim * col_dim, 1.0, tmp_mat, 1, ret, 1);

	// 解放
	free(high_mat_a);
	free(low_mat_a);
	free(high_mat_b);
	free(low_mat_b);
	//free(tmp_mat);
}

// matrix * vector
void accurate_dgemv(double ret[], double mat[], int row_dim, int col_dim, double vec[])
{
	int i, j, ij;
	double p, s, q, r, h;

	for(i = 0; i < row_dim; i++)
	{
		//ret[i] = 0;
		p = 0; s = 0;
		for(j = 0; j < col_dim; j++)
		{
			ij = i * col_dim + j;

			//ret[i] += mat[ij] * vec[j];
			h = two_prod(mat[ij], vec[j], &r);
			p = two_sum(p, h, &q);
			s += q + r;
		}
		ret[i] = p + s;
	}
}



// SplitMat_A (Lower triangular matrix)
// ltrimat = [ a_11 .................. ]
//           [ a_21 a_22 ............. ]
//           [ ....................... ]
//           [ a_m1 a_m2 ........ a_mn ]
void dsplit_ltrimat(double high_mat[], double low_mat[], double ltrimat[], int row_dim, int col_dim)
{
	int i, j;
	int num_digits = 53; // IEEE double prec.
	double *s;
	double mu, abs_aij, t_exp;

	// 初期化
	s = (double *)calloc(row_dim * col_dim, sizeof(double));

	// s := 0
	drepeat_constant_array(s, 0.0, row_dim * col_dim);

	// mu[i] = max_j |mat[i, j]|
	for(i = 0; i < row_dim; i++)
	{
		mu = fabs(ltrimat[i * col_dim + 0]);
		//for(j = 1; j < col_dim; j++)
		for(j = 1; j <= i; j++)
		{
			abs_aij = fabs(ltrimat[i * col_dim + j]);
			if(abs_aij > mu)
				mu = abs_aij;
		}

		// t_exp = ceil(log2(mu)) + ceil(s + log2(col_dim + 1) / 2)
		t_exp = ceil(DLOG2(mu)) + ceil(((double)num_digits + DLOG2((double)(col_dim + 1))/ 2.0));

		// s[i, j] = 2^t_exp
		//for(j = 0; j < col_dim; j++)
		for(j = 0; j <= i; j++)
			s[i * col_dim + j] = pow(2.0, t_exp);
	}

	// tmp_mat := mat + s
	cblas_dcopy(row_dim * col_dim, ltrimat, 1, high_mat, 1);
	cblas_daxpy(row_dim * col_dim, 1.0, s, 1, high_mat, 1);

	// high_mat := tmp_mat - s
	cblas_daxpy(row_dim * col_dim, -1.0, s, 1, high_mat, 1);

	// low_mat := mat - high_mat
	cblas_dcopy(row_dim * col_dim, ltrimat, 1, low_mat, 1);
	cblas_daxpy(row_dim * col_dim, -1.0, high_mat, 1, low_mat, 1);

	// 解放
	free(s);
}

// SplitMat_A (Upper triangular matrix)
// utrimat = [ a_11 a_12 ........ a_1n ]
//           [      a_22 ........ a_2n ]
//           [           ............. ]
//           [                    a_mn ]
void dsplit_utrimat(double high_mat[], double low_mat[], double utrimat[], int row_dim, int col_dim)
{
	int i, j;
	int num_digits = 53; // IEEE double prec.
	double *s;
	double mu, abs_aij, t_exp;

	// 初期化
	s = (double *)calloc(row_dim * col_dim, sizeof(double));

	// s := 0
	drepeat_constant_array(s, 0.0, row_dim * col_dim);

	// mu[i] = max_j |mat[i, j]|
	for(i = 0; i < row_dim; i++)
	{
		mu = fabs(utrimat[i * col_dim + i]);
		//for(j = 1; j < col_dim; j++)
		for(j = i + 1; j < col_dim; j++)
		{
			abs_aij = fabs(utrimat[i * col_dim + j]);
			if(abs_aij > mu)
				mu = abs_aij;
		}

		// t_exp = ceil(log2(mu)) + ceil(s + log2(col_dim + 1) / 2)
		t_exp = ceil(DLOG2(mu)) + ceil(((double)num_digits + DLOG2((double)(col_dim + 1)) / 2.0));

		// s[i, j] = 2^t_exp
		//for(j = 0; j < col_dim; j++)
		for(j = i; j < col_dim; j++)
			s[i * col_dim + j] = pow(2.0, t_exp);
	}

	// tmp_mat := mat + s
	cblas_dcopy(row_dim * col_dim, utrimat, 1, high_mat, 1);
	cblas_daxpy(row_dim * col_dim, 1.0, s, 1, high_mat, 1);

	// high_mat := tmp_mat - s
	cblas_daxpy(row_dim * col_dim, -1.0, s, 1, high_mat, 1);

	// low_mat := mat - high_mat
	cblas_dcopy(row_dim * col_dim, utrimat, 1, low_mat, 1);
	cblas_daxpy(row_dim * col_dim, -1.0, high_mat, 1, low_mat, 1);

	// 解放
	free(s);
}


// Accurate_TRMM
// ret[row_dim * col_dim] := mat_tri[row_dim * mid_dim] * mat_b[mid_dim * col_dim]
void accurate_dtrmm(double ret[], double mat_tri[], CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans, CBLAS_DIAG diag, double mat_b[], int row_dim, int col_dim)
{
	double *high_mat_a, *low_mat_a, *high_mat_b, *low_mat_b, *tmp_mat;

	// 初期化
	high_mat_a = (double *)calloc(row_dim * col_dim, sizeof(double));
	low_mat_a = (double *)calloc(row_dim * col_dim, sizeof(double));
	high_mat_b = (double *)calloc(row_dim * col_dim, sizeof(double));
	low_mat_b = (double *)calloc(row_dim * col_dim, sizeof(double));
	tmp_mat = (double *)calloc(row_dim * col_dim, sizeof(double));

	// ret, tmp_mat := 0
	//drepeat_constant_array(ret, 0.0, row_dim * col_dim);
	//drepeat_constant_array(tmp_mat, 0.0, row_dim * col_dim);

	// 分割
	if(uplo == CblasUpper)
		dsplit_utrimat(high_mat_a, low_mat_a, mat_tri, row_dim, col_dim);
	else if(uplo == CblasLower)
		dsplit_ltrimat(high_mat_a, low_mat_a, mat_tri, row_dim, col_dim);
	else
	{
		fprintf(stderr, "ERROR: accurate_dtrmm can not receive %d as uplo!\n", uplo);

		// 解放
		free(high_mat_a);
		free(low_mat_a);
		free(high_mat_b);
		free(low_mat_b);
		free(tmp_mat);

		return;
	}

	//dsplit_mat(high_mat_a, low_mat_a, mat_tri, row_dim, col_dim);

	dsplit_mat_t(high_mat_b, low_mat_b, mat_b, row_dim, col_dim);

	// ret += high_mat_a * low_mat_b + low_mat_a * mat_b
	// ret += high_mat_a * high_mat_b
	cblas_dcopy(row_dim * col_dim, low_mat_b, 1, ret, 1);
	cblas_dtrmm(CblasRowMajor, CblasLeft, uplo, trans, diag, row_dim, col_dim, 1.0, high_mat_a, row_dim, ret, col_dim);
	cblas_dcopy(row_dim * col_dim, mat_b, 1, tmp_mat, 1);
	cblas_dtrmm(CblasRowMajor, CblasLeft, uplo, trans, diag, row_dim, col_dim, 1.0, low_mat_a, row_dim, tmp_mat, col_dim);
	cblas_daxpy(row_dim * col_dim, 1.0, tmp_mat, 1, ret, 1);
	cblas_dtrmm(CblasRowMajor, CblasLeft, uplo, trans, diag, row_dim, col_dim, 1.0, high_mat_a, row_dim, high_mat_b, col_dim);
	cblas_daxpy(row_dim * col_dim, 1.0, high_mat_b, 1, ret, 1);

	// 解放
	free(high_mat_a);
	free(low_mat_a);
	free(high_mat_b);
	free(low_mat_b);
	free(tmp_mat);
}

// triangular matrix * vector
void accurate_dtrmv(double ret[], double tri_mat[], CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans, int row_dim, int col_dim, double vec[])
{
	int i, j, ij;
	double p, s, q, r, h;

	// Upper Triangular matrix
	if(uplo == CblasUpper)
	{
		if(trans == CblasNoTrans)
		{
			for(i = 0; i < row_dim; i++)
			{
				//ret[i] = 0;
				//p = 0; s = 0;
				//for(j = i; j < col_dim; j++)
				p = two_prod(tri_mat[i * col_dim + i], vec[i], &s);

				for(j = i + 1; j < col_dim; j++)
				{
					ij = i * col_dim + j;

					//ret[i] += tri_mat[ij] * vec[j];
					h = two_prod(tri_mat[ij], vec[j], &r);
					p = two_sum(p, h, &q);
					s += q + r;
				}
				ret[i] = p + s;
			}
		}
		else // Transposed -> Lower Triangular matrix
		{
			for(i = 0; i < col_dim; i++)
			{
				//ret[i] = 0.0;
				p = 0; s = 0;
				//for(j = 0; j <= i ; j++)
				p = two_prod(tri_mat[i], vec[0], &s);
				for(j = 1; j <= i ; j++)
				{
					//ij = i * row_dim + j;
					ij = i + j * col_dim;

					//ret[i] += tri_mat[ij] * vec[j];
					h = two_prod(tri_mat[ij], vec[j], &r);
					p = two_sum(p, h, &q);
					s += q + r;
				}
				ret[i] = p + s;
			}
		}
	}

	// Lower Triangular matrix
	if(uplo == CblasLower)
	{
		if(trans == CblasNoTrans)
		{
			for(i = 0; i < row_dim; i++)
			{
				//ret[i] = 0.0;
				//p = 0; s = 0;
				//for(j = 0; j <= i; j++)
				p = two_prod(tri_mat[i * col_dim], vec[0], &s);
				for(j = 1; j <= i; j++)
				{
					ij = i * col_dim + j;

					//ret[i] += tri_mat[ij] * vec[j];
					h = two_prod(tri_mat[ij], vec[j], &r);
					p = two_sum(p, h, &q);
					s += q + r;
				}
				ret[i] = p + s;
			}
		}
		else // Transposed -> Upper Triangular matrix
		{
			for(i = 0; i < col_dim; i++)
			{
				//ret[i] = 0.0;
				//p = 0; s = 0;
				//for(j = i; j < row_dim; j++)
				p = two_prod(tri_mat[i + i * col_dim], vec[i], &s);
				for(j = i + 1; j < row_dim; j++)
				{
					//ij = i * row_dim + j;
					ij = i + j * col_dim;

					//ret[i] += tri_mat[ij] * vec[j];
					h = two_prod(tri_mat[ij], vec[j], &r);
					p = two_sum(p, h, &q);
					s += q + r;
				}
				ret[i] = p + s;
			}
		}
	}
}

// check pivot
int check_pivot(int pivot[], int dim)
{
	int miss = 0;
	while(--dim >= 0)
	{
		if(pivot[dim] > (dim + 1))
		{
			//printf("pivot[%d] != %d\n", dim, pivot[dim]);
			miss++;
		}
	}

	return miss;
}

// pivotに従った行の入れ替え
/*
pivot[0] = 10 (1) 9 <-> 0	9				9	
pivot[1] = 7				1	(2) 6 <-> 1	6	
pivot[2] = 5				2				2	(3) 4 <-> 2
pivot[3] = 8				3				3	(4) 7 <-> 3
pivot[4] = 7				4				4	(5) 6 <-> 7
pivot[5] = 10				5				5
pivot[6] = 10				6				1
pivot[7] = 9				7				7
pivot[8] = 9				8				8
pivot[9] = 10				0				0
*/
void pivot_order_dmat(double mat[], int row_dim, int col_dim, int pivot[])
{
	int i, j, *row_order, tmp_i;
	double *tmp_mat;

	row_order = (int *)calloc(row_dim, sizeof(int));
	tmp_mat = (double *)calloc(row_dim * col_dim, sizeof(double));

	for(i = 0; i < row_dim; i++)
		row_order[i] = i;

	for(i = 0; i < row_dim; i++)
	{
		if(pivot[i] > (i + 1))
		{
			tmp_i = row_order[i];
			row_order[i] = row_order[pivot[i] - 1];
			row_order[pivot[i] - 1] = tmp_i;
		}
	}

	//for(i = row_dim - 1; i >= 0; i--)
	for(i = 0; i < row_dim; i++)
	{
		//printf("pivot[%d] = %d\n", i, pivot[i]);
		for(j = 0; j < col_dim; j++)
			tmp_mat[i * col_dim + j] = mat[row_order[i] * col_dim + j];
	}

	cblas_dcopy(row_dim * col_dim, tmp_mat, 1, mat, 1);

	free(row_order);
	free(tmp_mat);
}
void inverse_pivot_order_dmat(double mat[], int row_dim, int col_dim, int pivot[])
{
	int i, j, *row_order, tmp_i;
	double *tmp_mat;

	row_order = (int *)calloc(row_dim, sizeof(int));
	tmp_mat = (double *)calloc(row_dim * col_dim, sizeof(double));

	for(i = 0; i < row_dim; i++)
		row_order[i] = i;

	for(i = 0; i < row_dim; i++)
	{
		if(pivot[i] > (i + 1))
		{
			tmp_i = row_order[i];
			row_order[i] = row_order[pivot[i] - 1];
			row_order[pivot[i] - 1] = tmp_i;
		}
	}

	//for(i = row_dim - 1; i >= 0; i--)
	for(i = 0; i < row_dim; i++)
	{
		//printf("pivot[%d] = %d\n", i, pivot[i]);
		for(j = 0; j < col_dim; j++)
			tmp_mat[row_order[i] * col_dim + j] = mat[i * col_dim + j];
			//tmp_mat[i * col_dim + j] = mat[row_order[i] * col_dim + j];
	}

	cblas_dcopy(row_dim * col_dim, tmp_mat, 1, mat, 1);

	free(row_order);
	free(tmp_mat);
}

// transpose
// ret := org_mat^T
void transpose_dmatrix(double ret[], double org_mat[], int row_dim, int col_dim)
{
	int i, j;

	for(i = 0; i < row_dim; i++)
		for(j = 0; j < col_dim; j++)
			ret[j * row_dim + i] = org_mat[i * col_dim + j];

}


// print_vec
void print_vec(double vec[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		printf("%5d, %25.17e\n", i, vec[i]);
}
// print_vec
void print_vec2(double vec[], double vec2[], int dim)
{
	int i;

	for(i = 0; i < dim; i++)
		printf("%5d, %25.17e, %25.17e\n", i, vec[i], vec2[i]);
}

// Iterative ref.
// ret[dim] := A[dim * dim]^(-1) * b[dim]
int d_iterative_refinement_amm(double ret[], double mat_a[], double vec_b[], int dim, int max_itimes, double rtol, double atol)
{
	lapack_int info, itimes, i;
	lapack_int *pivot;
	int acculate_flag = 0;
	double *vec_z, *vec_r, *vec_x, *vec_d, *mat_lu, *mat_c, *mat_at;
	double norm_z, old_norm_z, norm_x, alpha = 0.1;
	double cond_a, cond_c, norm_mat_a, norm_mat_c;

	// ピボット・ベクトル初期化
	pivot = (lapack_int *)calloc(dim, sizeof(lapack_int));
	mat_lu = (double *)calloc(dim * dim, sizeof(double));
	mat_at = (double *)calloc(dim * dim, sizeof(double));
	mat_c = (double *)calloc(dim * dim, sizeof(double));
	vec_z = (double *)calloc(dim, sizeof(double));
	vec_d = (double *)calloc(dim, sizeof(double));
	vec_r = (double *)calloc(dim, sizeof(double));
	vec_x = (double *)calloc(dim, sizeof(double));

	// mat_lu := mat_a
	cblas_dcopy(dim * dim, mat_a, 1, mat_lu, 1);

	// LU分解: P * A^T := L * U
	// mat_lu := mat_a^T
	transpose_dmatrix(mat_at, mat_a, dim, dim);
	cblas_dcopy(dim * dim, mat_at, 1, mat_lu, 1);
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_lu, dim, pivot);
	//printf("DGETRF info = %d\n", info);

	printf("pivot has been changed!(%d)\n", check_pivot(pivot, dim));

	// x0 := 0
	// r0 := b
	drepeat_constant_array(vec_x, 0.0, dim);

	// P^(-1) * r
	//inverse_pivot_order_dmat(mat_a, dim, dim, pivot);
	//inverse_pivot_order_dmat(vec_b, dim, 1, pivot);

	norm_x = cblas_dnrm2(dim, vec_x, 1);
	cblas_dcopy(dim, vec_b, 1, vec_r, 1);
	cblas_dcopy(dim, vec_r, 1, vec_z, 1);

	old_norm_z = cblas_dnrm2(dim, vec_z, 1);

	// 条件数
	norm_mat_a = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1', dim, dim, mat_a, dim);
	LAPACKE_dgecon(LAPACK_ROW_MAJOR, '1', dim, mat_lu, dim, norm_mat_a, &cond_a);
	cond_a = 1.0 / cond_a;

	// L^(-1), U^(-1)
	info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'L', 'U', dim, mat_lu, dim);
	info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'U', 'N', dim, mat_lu, dim);
	//printf("DTRTRI info = %d\n", info);

	// ループ
	for(itimes = 0; itimes < max_itimes; itimes++)
	{
		// mat_a * z = rを解く
		// 前進，後退代入
		//info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', dim, 1, mat_lu, dim, pivot, vec_z, 1);
		//printf("DGETRS info = %d\n", info);

		// z := P^(-1) * L^(-T) * U^(-T) * r
		cblas_dtrmv(CblasRowMajor, CblasUpper, CblasTrans, CblasNonUnit, dim, mat_lu, dim, vec_z, 1);
		cblas_dtrmv(CblasRowMajor, CblasLower, CblasTrans, CblasUnit, dim, mat_lu, dim, vec_z, 1);
		inverse_pivot_order_dmat(vec_z, dim, 1, pivot);

		// x += z
		cblas_daxpy(dim, 1.0, vec_z, 1, vec_x, 1);

		print_vec2(vec_x, vec_z, dim);

		// 収束判定(1)
		// ||(x + z) - x|| = ||z|| <= rtol * ||x|| + atol ?
		norm_z = cblas_dnrm2(dim, vec_z, 1);
		printf("%5d, %25.17e\n", itimes, norm_z);
		if(norm_z <= rtol * norm_x + atol)
			break;


		// 収束判定(2)
		// ||z|| >= alpha * ||old_z||
		if(norm_z >= alpha * old_norm_z)
		{
			if(acculate_flag != 0)
				break;

			// pivotに従った行の入れ替え
			//inverse_pivot_order_dmat(mat_a, dim, dim, pivot);
			//pivot_order_dmat(mat_at, dim, dim, pivot);
			//pivot_order_dmat(mat_a, dim, dim, pivot);
			//pivot_order_dmat(mat_lu, dim, dim, pivot);
			//pivot_order_dmat(vec_b, dim, 1, pivot);

			//transpose_dmatrix(mat_a, mat_at, dim, dim);

			// U^(-1)
			//info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'U', 'N', dim, mat_lu, dim);
			//printf("DTRTRI info = %d\n", info);
			
			// C := U^(-1)^T * A
			//cblas_dcopy(dim * dim, mat_a, 1, mat_c, 1);
			//cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, dim, dim, 1.0, mat_lu, dim, mat_c, dim);
			accurate_dtrmm(mat_c, mat_lu, CblasUpper, CblasTrans, CblasNonUnit, mat_a, dim, dim);
			//accurate_dtrmm(mat_c, mat_lu, CblasUpper, CblasNoTrans, CblasNonUnit, mat_a, dim, dim);
			//accurate_dgemm(mat_c, mat_lu, mat_a, dim, dim, dim);

			// ||C||
			norm_mat_c = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1', dim, dim, mat_c, dim);

			// d := U^(-1)^T * b
			accurate_dtrmv(vec_d, mat_lu, CblasUpper, CblasTrans, dim, dim, vec_b);
			//accurate_dtrmv(vec_d, mat_lu, CblasUpper, CblasNoTrans, dim, dim, vec_b);

			// A := C
			cblas_dcopy(dim * dim, mat_c, 1, mat_a, 1);

			// mat_lu := mat_a^T
			transpose_dmatrix(mat_lu, mat_c, dim, dim);

			// LU := C
			//cblas_dcopy(dim * dim, mat_c, 1, mat_lu, 1);
			info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_lu, dim, pivot);

			// 条件数
			LAPACKE_dgecon(LAPACK_ROW_MAJOR, '1', dim, mat_lu, dim, norm_mat_c, &cond_c);
			cond_c = 1.0 / cond_c;

			// L^(-1), U^(-1)
			info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'L', 'U', dim, mat_lu, dim);
			info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'U', 'N', dim, mat_lu, dim);

			// b := d
			cblas_dcopy(dim, vec_d, 1, vec_b, 1);

			//check_pivot(pivot, dim);
			printf("A := C! cond1(A) = %15.7e -> cond1(C) = %15.7e\n", cond_a, cond_c);
					
			// A has been changed!
			acculate_flag = 1;
		}
	
		old_norm_z = norm_z;

		// ||x||
		norm_x = cblas_dnrm2(dim, vec_x, 1);

		// r := vec_b - mat_a * x
		cblas_dcopy(dim, vec_b, 1, vec_r, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, -1.0, mat_a, dim, vec_x, 1, 1.0, vec_r, 1);
		//accurate_dgemv(vec_r, mat_a, dim, dim, vec_x);
		//cblas_daxpy(dim, -1.0, vec_b, 1, vec_r, 1);
		//cblas_dscal(dim, -1.0, vec_r, 1);

		// z := r
		cblas_dcopy(dim, vec_r, 1, vec_z, 1);
	}

	// ret := x
	cblas_dcopy(dim, vec_x, 1, ret, 1);

	// ピボット・ベクトル消去
	free(pivot);
	free(mat_lu);
	free(mat_at);
	free(vec_z);
	free(vec_d);
	free(vec_r);
	free(vec_x);

	return itimes;
}

// Iterative ref.
// ret[dim] := A[dim * dim]^(-1) * b[dim]
int d_iterative_refinement_amm2(double ret[], double mat_a[], double vec_b[], int dim, int max_itimes, double rtol, double atol)
{
	lapack_int info, itimes, i;
	lapack_int *pivot;
	int acculate_flag = 0;
	double *vec_z, *vec_r, *vec_x, *vec_d, *mat_lu, *mat_c, *mat_at, *org_mat_a, *org_vec_b;
	double norm_z, old_norm_z, norm_x, alpha = 0.1;
	double cond_a, cond_c, norm_mat_a, norm_mat_c;

	// ピボット・ベクトル初期化
	pivot = (lapack_int *)calloc(dim, sizeof(lapack_int));
	mat_lu = (double *)calloc(dim * dim, sizeof(double));
	org_mat_a = (double *)calloc(dim * dim, sizeof(double));
	mat_at = (double *)calloc(dim * dim, sizeof(double));
	mat_c = (double *)calloc(dim * dim, sizeof(double));
	vec_z = (double *)calloc(dim, sizeof(double));
	org_vec_b = (double *)calloc(dim, sizeof(double));
	vec_d = (double *)calloc(dim, sizeof(double));
	vec_r = (double *)calloc(dim, sizeof(double));
	vec_x = (double *)calloc(dim, sizeof(double));

	// mat_lu := mat_a
	cblas_dcopy(dim * dim, mat_a, 1, mat_lu, 1);
	cblas_dcopy(dim * dim, mat_a, 1, org_mat_a, 1);

	// LU分解: P * A^T := L * U
	// mat_lu := mat_a^T
	transpose_dmatrix(mat_at, mat_a, dim, dim);
	cblas_dcopy(dim * dim, mat_at, 1, mat_lu, 1);
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_lu, dim, pivot);
	//printf("DGETRF info = %d\n", info);

	printf("pivot has been changed!(%d)\n", check_pivot(pivot, dim));

	// 条件数
	norm_mat_a = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1', dim, dim, mat_a, dim);
	LAPACKE_dgecon(LAPACK_ROW_MAJOR, '1', dim, mat_lu, dim, norm_mat_a, &cond_a);
	cond_a = 1.0 / cond_a;

	// L^(-1), U^(-1)
	info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'L', 'U', dim, mat_lu, dim);
	info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'U', 'N', dim, mat_lu, dim);
	//printf("DTRTRI info = %d\n", info);

	if(cond_a >= 1.0e+5)
	{
		// pivotに従った行の入れ替え
		//inverse_pivot_order_dmat(mat_a, dim, dim, pivot);
		//pivot_order_dmat(mat_at, dim, dim, pivot);
		//pivot_order_dmat(mat_a, dim, dim, pivot);
		//pivot_order_dmat(mat_lu, dim, dim, pivot);
		//pivot_order_dmat(vec_b, dim, 1, pivot);

		//transpose_dmatrix(mat_a, mat_at, dim, dim);

		// U^(-1)
		//info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'U', 'N', dim, mat_lu, dim);
		//printf("DTRTRI info = %d\n", info);
		
		// C := U^(-1)^T * A
		//cblas_dcopy(dim * dim, mat_a, 1, mat_c, 1);
		//cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, dim, dim, 1.0, mat_lu, dim, mat_c, dim);
		accurate_dtrmm(mat_c, mat_lu, CblasUpper, CblasTrans, CblasNonUnit, mat_a, dim, dim);
		//accurate_dtrmm(mat_c, mat_lu, CblasUpper, CblasNoTrans, CblasNonUnit, mat_a, dim, dim);
		//accurate_dgemm(mat_c, mat_lu, mat_a, dim, dim, dim);

		// ||C||
		norm_mat_c = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1', dim, dim, mat_c, dim);

		// d := U^(-1)^T * b
		accurate_dtrmv(vec_d, mat_lu, CblasUpper, CblasTrans, dim, dim, vec_b);
		//accurate_dtrmv(vec_d, mat_lu, CblasUpper, CblasNoTrans, dim, dim, vec_b);

		// A := C
		cblas_dcopy(dim * dim, mat_c, 1, mat_a, 1);

		// mat_lu := mat_a^T
		transpose_dmatrix(mat_lu, mat_c, dim, dim);

		// LU := C
		//cblas_dcopy(dim * dim, mat_c, 1, mat_lu, 1);
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_lu, dim, pivot);

		// 条件数
		LAPACKE_dgecon(LAPACK_ROW_MAJOR, '1', dim, mat_lu, dim, norm_mat_c, &cond_c);
		cond_c = 1.0 / cond_c;

		// L^(-1), U^(-1)
		info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'L', 'U', dim, mat_lu, dim);
		info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'U', 'N', dim, mat_lu, dim);

		// b := d
		cblas_dcopy(dim, vec_d, 1, vec_b, 1);

		//check_pivot(pivot, dim);
		printf("A := C! cond1(A) = %15.7e -> cond1(C) = %15.7e\n", cond_a, cond_c);
				
		// A has been changed!
		acculate_flag = 1;
	}

	// x0 := 0
	// r0 := b
	drepeat_constant_array(vec_x, 0.0, dim);

	// P^(-1) * r
	//inverse_pivot_order_dmat(mat_a, dim, dim, pivot);
	//inverse_pivot_order_dmat(vec_b, dim, 1, pivot);

	norm_x = cblas_dnrm2(dim, vec_x, 1);
	cblas_dcopy(dim, vec_b, 1, vec_r, 1);
	cblas_dcopy(dim, vec_r, 1, vec_z, 1);

	old_norm_z = cblas_dnrm2(dim, vec_z, 1);

	// ループ
	for(itimes = 0; itimes < max_itimes; itimes++)
	{
		// mat_a * z = rを解く
		// 前進，後退代入
		//info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', dim, 1, mat_lu, dim, pivot, vec_z, 1);
		//printf("DGETRS info = %d\n", info);

		// z := P^(-1) * L^(-T) * U^(-T) * r
		cblas_dtrmv(CblasRowMajor, CblasUpper, CblasTrans, CblasNonUnit, dim, mat_lu, dim, vec_z, 1);
		cblas_dtrmv(CblasRowMajor, CblasLower, CblasTrans, CblasUnit, dim, mat_lu, dim, vec_z, 1);
		inverse_pivot_order_dmat(vec_z, dim, 1, pivot);

		// x += z
		cblas_daxpy(dim, 1.0, vec_z, 1, vec_x, 1);

		print_vec2(vec_x, vec_z, dim);

		// 収束判定(1)
		// ||(x + z) - x|| = ||z|| <= rtol * ||x|| + atol ?
		norm_z = cblas_dnrm2(dim, vec_z, 1);
		printf("%5d, %25.17e\n", itimes, norm_z);
		if(norm_z <= rtol * norm_x + atol)
			break;


		// 収束判定(2)
		// ||z|| >= alpha * ||old_z||
		//if(norm_z >= alpha * old_norm_z)
		//	break;
	
		old_norm_z = norm_z;

		// ||x||
		norm_x = cblas_dnrm2(dim, vec_x, 1);

		// r := U * ((vec_b - mat_a * x)
		//cblas_dcopy(dim, org_vec_b, 1, vec_r, 1);
		//cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, -1.0, org_mat_a, dim, vec_x, 1, 1.0, vec_r, 1);
		//cblas_dtrmv(CblasRowMajor, CblasUpper, CblasTrans, CblasNonUnit, dim, mat_lu, dim, vec_r, 1);

		accurate_dgemv(vec_r, mat_a, dim, dim, vec_x);
		cblas_daxpy(dim, -1.0, vec_b, 1, vec_r, 1);
		cblas_dscal(dim, -1.0, vec_r, 1);

		// z := r
		cblas_dcopy(dim, vec_r, 1, vec_z, 1);
	}

	// ret := x
	cblas_dcopy(dim, vec_x, 1, ret, 1);

	// ピボット・ベクトル消去
	free(pivot);
	free(mat_lu);
	free(mat_at);
	free(org_mat_a);
	free(vec_z);
	free(org_vec_b);
	free(vec_d);
	free(vec_r);
	free(vec_x);

	return itimes;
}

// MPFR
#include "mpfr.h"
#include "linear_c.h"

#ifdef __MPFR_H
// ret := source
void set_array_d2mpfr(mpfr_t ret[], double source[], int dim)
{
	while(--dim >= 0)
		mpfr_set_d(ret[dim], source[dim], _tk_default_rmode);
}
void set_array_mpfr2d(double ret[], mpfr_t source[], int dim)
{
	while(--dim >= 0)
		ret[dim] = mpfr_get_d(source[dim], _tk_default_rmode);
}

void precise_dgemv(double ret[], double mat[], double x[], int dim)
{
	unsigned long prec = 256;
	mpfr_t *mpfr_mat, *mpfr_ret, *mpfr_x;

	mpfr_mat = (mpfr_t *)calloc(dim * dim, sizeof(mpfr_t));
	mpfr_ret = (mpfr_t *)calloc(dim, sizeof(mpfr_t));
	mpfr_x = (mpfr_t *)calloc(dim, sizeof(mpfr_t));

	mpfr_init2_array(mpfr_mat, dim * dim, prec);
	mpfr_init2_array(mpfr_ret, dim, prec);
	mpfr_init2_array(mpfr_x, dim, prec);

	set_array_d2mpfr(mpfr_mat, mat, dim * dim);
	set_array_d2mpfr(mpfr_x, x, dim);

	mpfr_mymv(mpfr_ret, mpfr_mat, mpfr_x, dim);

	set_array_mpfr2d(ret, mpfr_ret, dim);

	mpfr_clear_array(mpfr_mat, dim * dim);
	mpfr_clear_array(mpfr_ret, dim);
	mpfr_clear_array(mpfr_x, dim);
	free(mpfr_mat);
	free(mpfr_ret);
	free(mpfr_x);
}
#endif // __MPFR_H

// print_dmatrix
void print_dmatrix(double mat[], int pivot[], int row_dim, int col_dim)
{
	int i, j;

	for(i = 0; i < row_dim; i++)
	{
		if(pivot != NULL)
			printf("%3d %3d: ", i, pivot[i]);
		else
			printf("%3d: ", i);
		
		for(j = 0; j < col_dim; j++)
			printf("%12.5e ", mat[i * col_dim + j]);
		printf("\n");
	}
}

// setI
void setI_dmatrix(double mat[], int row_dim, int col_dim)
{
	int i, j;

	for(i = 0; i < row_dim; i++)
	{
		for(j = 0; j < col_dim; j++)
			mat[i * col_dim + j] = 0.0;

		if(i < col_dim)
			mat[i * col_dim + i] = 1.0;
	}
}

int main()
{
	lapack_int i, j, dim, itimes;
	lapack_int inc_vec_x, inc_vec_b;
	lapack_int *pivot, info;

	double *mat_a, *vec_b, *vec_x, *vec_ret, *mat_at;
	double *org_mat_a, *org_vec_b, *tmp_mat;
	double alpha, beta, norm_a;
	double running_time;

	double test5x5_a[] = {
	   0.123456428858216,   0.115598209785209,  -0.277660682277346,   0.084040431208059,  -0.322377486163418, \
	   0.219541368533906,   0.203546499742431,  -0.491390560544294,   0.148720453439147,  -0.572721901709023, \
	   0.033255324299145,   0.031189739335185,  -0.074852124886012,   0.022657639724993,  -0.086854172143882, \
	  -0.070510178088673,  -0.068196155097055,   0.161116663463259,  -0.048800625091667,   0.184745468471442, \
	  -0.031350857196837,  -0.030368613494423,   0.071693731663506,  -0.021712612295286,   0.082153206459781 \
	};
	double test5x5_b[] = {  -0.276943098589280,  -0.492304140537833, -0.074603593670572, 0.158355173657306, 0.070414855136741 };
	
	// 次元数入力
	printf("Dim = "); scanf("%d", &dim);

	if(dim <= 0)
	{
		printf("Illigal dimenstion! (dim = %d)\n", dim);
		return EXIT_FAILURE;
	}

	// 変数初期化
	mat_a = (double *)calloc(dim * dim, sizeof(double));
	mat_at = (double *)calloc(dim * dim, sizeof(double));
	org_mat_a = (double *)calloc(dim * dim, sizeof(double));
	tmp_mat = (double *)calloc(dim * dim, sizeof(double));
	vec_x = (double *)calloc(dim, sizeof(double));
	vec_b = (double *)calloc(dim, sizeof(double));
	org_vec_b = (double *)calloc(dim, sizeof(double));
	vec_ret = (double *)calloc(dim, sizeof(double));
	pivot = (lapack_int *)calloc(dim, sizeof(lapack_int));

	// mat_aとvec_xに値入力
	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
		{
			//mat_a[i * dim + j] = (double)(dim - ((i > j) ? i : j));
			//mat_a[i * dim + j] = (double)rand();
			//mat_a[i * dim + j] = (double)rand() / (double)RAND_MAX;
			//mat_a[i * dim + j] = 1.0 / (double)(i + j + 1);
			//mat_a[i * dim + j] = test5x5_a[i * dim + j];
			//if((i + j + 1) % 2 != 0)
			//	mat_a[i * dim + j] = -mat_a[i * dim + j];

			mat_a[i * dim + j] = 1.0 / (double)(i + j + 1);
		}
		//mat_a[i * dim + i] += 2.0;
		//vec_x[i] = 1.0 / (double)(i + 1);
		//vec_x[i] = (double)(i + 1);
		vec_x[i] = 1.0;
	}

	// size(vec_x) == size(vec_b)
	inc_vec_x = inc_vec_b = 1;

	// vec_b := 1.0 * mat_a * vec_x + 0.0 * vec_b
	alpha = 1.0;
	beta = 0.0;
	//cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, alpha, mat_a, dim, vec_x, inc_vec_x, beta, vec_b, inc_vec_b);
	precise_dgemv(vec_b, mat_a, vec_x, dim);

	// org_mat_a := mat_a
	// org_vec_b := vec_b
	cblas_dcopy(dim * dim, mat_a, 1, org_mat_a, 1);
	cblas_dcopy(dim, vec_b, 1, org_vec_b, 1);

	// 出力
/*	for(i = 0; i < dim; i++)
	{
		printf("[");
		for(j = 0; j < dim; j++)
			printf("%10.3f ", mat_a[i * dim + j]);
		printf("]  %10.3f = %10.3f\n", vec_x[i], vec_b[i]);
	}
*/
	// 直接法
	printf("Direct Methods:%d x %d\n", dim, dim);
	cblas_dcopy(dim, vec_b, 1, vec_ret, 1);
	itimes = 0;

	//printf("mat_a: \n"); print_dmatrix(mat_a, NULL, dim, dim);
	
	//LAPACKE_dgesv(LAPACK_ROW_MAJOR, dim, 1, mat_a, dim, pivot, vec_ret, 1);
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_a, dim, pivot);
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', dim, 1, mat_a, dim, pivot, vec_ret, 1);

	//print_dmatrix(mat_a, pivot, dim, dim);

	// L * U * I
	setI_dmatrix(tmp_mat, dim, dim);
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, dim, dim, 1.0, mat_a, dim, tmp_mat, dim);
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, dim, dim, 1.0, mat_a, dim, tmp_mat, dim);

	/*
	printf("pivot: ");
	for(i = 0; i < dim; i++)
		printf("%d ", pivot[i]);
	printf("\n");
	*/
	
	//printf("tmp_mat: \n"); print_dmatrix(tmp_mat, pivot, dim, dim);

	// P * A - L * U
	cblas_dcopy(dim * dim, org_mat_a, 1, mat_a, 1);
	pivot_order_dmat(mat_a, dim, dim, pivot);

	//printf("mat_a: \n"); print_dmatrix(mat_a, pivot, dim, dim);

	cblas_daxpy(dim * dim, -1.0, tmp_mat, 1, mat_a, 1);

	norm_a = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1', dim, dim, org_mat_a, dim);
	norm_a = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1', dim, dim, mat_a, dim) / norm_a;
	printf("|| P * A - L * U ||_1 / ||A|| = %15.7e\n", norm_a);
	//print_dmatrix(mat_a, pivot, dim, dim);

	// diff
	//printf("x - calculated x = \n");
	for(i = 0; i < dim; i++)
	{
		vec_b[i] = fabs((vec_x[i] - vec_ret[i]) / vec_x[i]);
		//printf("%3d: ", i);
		//printf("%10.2e ", fabs((vec_x[i] - vec_ret[i]) / vec_x[i]));
		//printf("\n");
	}
	printf("||relerr(calculated x)||_2 = %25.17e\n", cblas_dnrm2(dim, vec_b, 1));

	fpu_fix_start(NULL);

	// 反復改良法
	printf("Iterative_ref Methods:%d x %d\n", dim, dim);
	cblas_dcopy(dim * dim, org_mat_a, 1, mat_a, 1);
	cblas_dcopy(dim, org_vec_b, 1, vec_b, 1);

	//itimes = d_iterative_refinement(vec_ret, mat_a, vec_b, dim, 10, 1.0e-15, 0.0);
	//itimes = d_iterative_refinement_amm(vec_ret, mat_a, vec_b, dim, 10, 1.0e-15, 0.0);
	//itimes = d_iterative_refinement_amm(vec_ret, mat_a, vec_b, dim, 10, 0.0, 0.0);
	itimes = d_iterative_refinement_amm2(vec_ret, mat_a, vec_b, dim, 10, 0.0, 0.0);

	// print
	printf("iteration times = %d\n", itimes);
	printf("calculated x = \n");
	for(i = 0; i < dim; i++)
	{
		//printf("%3d -> %3d: ", i, pivot[i]);
		printf("%25.17e ", vec_ret[i]);
		printf("\n");
	}

	// diff
	//printf("x - calculated x = \n");
	for(i = 0; i < dim; i++)
	{
		vec_b[i] = fabs((vec_x[i] - vec_ret[i]) / vec_x[i]);
		//printf("%3d: ", i);
		//printf("%10.2e ", fabs((vec_x[i] - vec_ret[i]) / vec_x[i]));
		//printf("\n");
	}
	printf("||relerr(calculated x)||_2 = %25.17e\n", cblas_dnrm2(dim, vec_b, 1));
	
	// 領域開放
	free(mat_a);
	free(mat_at);
	free(org_mat_a);
	free(tmp_mat);
	free(vec_x);
	free(vec_b);
	free(vec_ret);
	free(pivot);

	return EXIT_SUCCESS;
}
