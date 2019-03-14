/*************************************************/
/* linear_eq_dgetrf.c : �A���ꎟ����������
/* icc linear_eq.c -L/usr/local/lib -llapacke -llapack -lcblas -lrefblas -L/opt/intel/lib/intel64 -lifcore
/*************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "lapacke/lapacke.h"
#include "lapacke.h"
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

	// �s�{�b�g�E�x�N�g��������
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

	// LU����
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_lu, dim, pivot);
	//printf("DGETRF info = %d\n", info);

	// ���[�v
	for(itimes = 0; itimes < max_itimes; itimes++)
	{
		// mat_a * z = r������
		// �O�i�C��ޑ��
		info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', dim, 1, mat_lu, dim, pivot, vec_z, 1);
		//printf("DGETRS info = %d\n", info);

		// x += z
		cblas_daxpy(dim, 1.0, vec_z, 1, vec_x, 1);

		// ��������(1)
		// ||(x + z) - x|| = ||z|| <= rtol * ||x|| + atol ?
		norm_z = cblas_dnrm2(dim, vec_z, 1);
		if(norm_z <= rtol * norm_x + atol)
			break;

		// ��������(2)
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

	// �s�{�b�g�E�x�N�g������
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

	// ������
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
		t_exp = ceil(DLOG2(mu)) + ceil(((double)num_digits + DLOG2((double)(col_dim + 1))) / 2.0);

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

	// ���
	free(s);
}

// SplitMat_B
void dsplit_mat_t(double high_mat[], double low_mat[], double mat[], int row_dim, int col_dim)
{
	int i, j;
	int num_digits = 53; // IEEE double prec.
	double *s;
	double mu, abs_aij, t_exp;

	// ������
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
		t_exp = ceil(DLOG2(mu)) + ceil(((double)num_digits + DLOG2((double)(row_dim + 1))) / 2.0);

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

	// �o��
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
	// ���
	free(s);
}


// Accurate_MM
// ret[row_dim * col_dim] := mat_a[row_dim * mid_dim] * mat_b[mid_dim * col_dim]
void accurate_dgemm(double ret[], double mat_a[], double mat_b[], int row_dim, int col_dim, int mid_dim)
{
	double *high_mat_a, *low_mat_a, *high_mat_b, *low_mat_b, *tmp_mat;

	// ������
	high_mat_a = (double *)calloc(row_dim * mid_dim, sizeof(double));
	low_mat_a = (double *)calloc(row_dim * mid_dim, sizeof(double));
	high_mat_b = (double *)calloc(mid_dim * col_dim, sizeof(double));
	low_mat_b = (double *)calloc(mid_dim * col_dim, sizeof(double));
	//tmp_mat = (double *)calloc(row_dim * col_dim, sizeof(double));

	// ret, tmp_mat := 0
	drepeat_constant_array(ret, 0.0, row_dim * col_dim);
	drepeat_constant_array(tmp_mat, 0.0, row_dim * col_dim);

	// ����
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

	// ���
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

	// ������
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
		t_exp = ceil(DLOG2(mu)) + ceil(((double)num_digits + DLOG2((double)(col_dim + 1))) / 2.0);

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

	// ���
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

	// ������
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
		t_exp = ceil(DLOG2(mu)) + ceil(((double)num_digits + DLOG2((double)(col_dim + 1))) / 2.0);

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

	// ���
	free(s);
}


// Accurate_TRMM
// ret[row_dim * col_dim] := mat_tri[row_dim * mid_dim] * mat_b[mid_dim * col_dim]
void accurate_dtrmm(double ret[], double mat_tri[], CBLAS_UPLO uplo, CBLAS_TRANSPOSE trans, CBLAS_DIAG diag, double mat_b[], int row_dim, int col_dim)
{
	double *high_mat_a, *low_mat_a, *high_mat_b, *low_mat_b, *tmp_mat;

	// ������
	high_mat_a = (double *)calloc(row_dim * col_dim, sizeof(double));
	low_mat_a = (double *)calloc(row_dim * col_dim, sizeof(double));
	high_mat_b = (double *)calloc(row_dim * col_dim, sizeof(double));
	low_mat_b = (double *)calloc(row_dim * col_dim, sizeof(double));
	//tmp_mat = (double *)calloc(row_dim * col_dim, sizeof(double));

	// ret, tmp_mat := 0
	//drepeat_constant_array(ret, 0.0, row_dim * col_dim);
	//drepeat_constant_array(tmp_mat, 0.0, row_dim * col_dim);

	// ����
	if(uplo == CblasUpper)
		dsplit_utrimat(high_mat_a, low_mat_a, mat_tri, row_dim, col_dim);
	else if(uplo == CblasLower)
		dsplit_ltrimat(high_mat_a, low_mat_a, mat_tri, row_dim, col_dim);
	else
	{
		fprintf(stderr, "ERROR: accurate_dtrmm can not receive %d as uplo!\n", uplo);

		// ���
		free(high_mat_a);
		free(low_mat_a);
		free(high_mat_b);
		free(low_mat_b);
		free(tmp_mat);

		return;
	}

	dsplit_mat_t(high_mat_b, low_mat_b, mat_b, col_dim, col_dim);

	// ret += high_mat_a * low_mat_b + low_mat_a * mat_b
	// ret += high_mat_a * high_mat_b
	cblas_dtrmm(CblasRowMajor, CblasLeft, uplo, trans, diag, row_dim, col_dim, 1.0, high_mat_a, row_dim, low_mat_b, col_dim);
	cblas_dcopy(row_dim * col_dim, mat_b, 1, ret, 1);
	cblas_dtrmm(CblasRowMajor, CblasLeft, uplo, trans, diag, row_dim, col_dim, 1.0, low_mat_a, row_dim, ret, col_dim);
	cblas_daxpy(row_dim * col_dim, 1.0, low_mat_b, 1, ret, 1);
	cblas_dtrmm(CblasRowMajor, CblasLeft, uplo, trans, diag, row_dim, col_dim, 1.0, high_mat_a, row_dim, high_mat_b, col_dim);
	cblas_daxpy(row_dim * col_dim, 1.0, high_mat_b, 1, ret, 1);

	// ���
	free(high_mat_a);
	free(low_mat_a);
	free(high_mat_b);
	free(low_mat_b);
	//free(tmp_mat);
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
				p = 0; s = 0;
				for(j = i; j < col_dim; j++)
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
				for(j = 0; j <= i ; j++)
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
				p = 0; s = 0;
				for(j = 0; j <= i; j++)
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
				p = 0; s = 0;
				for(j = i; j < row_dim; j++)
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

// Iterative ref.
// ret[dim] := A[dim * dim]^(-1) * b[dim]
int d_iterative_refinement_amm(double ret[], double mat_a[], double vec_b[], int dim, int max_itimes, double rtol, double atol)
{
	lapack_int info, itimes, i;
	lapack_int *pivot;
	int acculate_flag = 0;
	double *vec_z, *vec_r, *vec_x, *vec_d, *mat_lu, *mat_c;
	double norm_z, old_norm_z, norm_x, alpha = 0.1;
	double cond_a, cond_c, norm_mat_a, norm_mat_c;

	// �s�{�b�g�E�x�N�g��������
	pivot = (lapack_int *)calloc(dim, sizeof(lapack_int));
	mat_lu = (double *)calloc(dim * dim, sizeof(double));
	mat_c = (double *)calloc(dim * dim, sizeof(double));
	vec_z = (double *)calloc(dim, sizeof(double));
	vec_d = (double *)calloc(dim, sizeof(double));
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

	// LU����
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_lu, dim, pivot);
	//printf("DGETRF info = %d\n", info);

	// ������
	norm_mat_a = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1', dim, dim, mat_a, dim);
	LAPACKE_dgecon(LAPACK_ROW_MAJOR, '1', dim, mat_lu, dim, norm_mat_a, &cond_a);
	cond_a = 1.0 / cond_a;

	// ���[�v
	for(itimes = 0; itimes < max_itimes; itimes++)
	{
		// mat_a * z = r������
		// �O�i�C��ޑ��
		info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', dim, 1, mat_lu, dim, pivot, vec_z, 1);
		//printf("DGETRS info = %d\n", info);

		// x += z
		cblas_daxpy(dim, 1.0, vec_z, 1, vec_x, 1);

		// ��������(1)
		// ||(x + z) - x|| = ||z|| <= rtol * ||x|| + atol ?
		norm_z = cblas_dnrm2(dim, vec_z, 1);
		printf("%5d, %25.17e\n", itimes, norm_z);
		if(norm_z <= rtol * norm_x + atol)
			break;


		// ��������(2)
		// ||z|| >= alpha * ||old_z||
		if(norm_z >= alpha * old_norm_z)
		{
			if(acculate_flag != 0)
				break;

			// U^(-1)
			LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'U', 'N', dim, mat_lu, dim);
			
			// C := U^(-1)^T * A
			accurate_dtrmm(mat_c, mat_lu, CblasUpper, CblasTrans, CblasNonUnit, mat_a, dim, dim);

			// ||C||
			norm_mat_c = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1', dim, dim, mat_c, dim);

			// d := U^(-1)^T * b
			accurate_dtrmv(vec_d, mat_lu, CblasUpper, CblasTrans, dim, dim, vec_b);

			// A := C
			cblas_dcopy(dim * dim, mat_c, 1, mat_a, 1);

			// LU := C
			cblas_dcopy(dim * dim, mat_c, 1, mat_lu, 1);
			info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, mat_lu, dim, pivot);

			// ������
			LAPACKE_dgecon(LAPACK_ROW_MAJOR, '1', dim, mat_lu, dim, norm_mat_a, &cond_c);
			cond_c = 1.0 / cond_c;

			// b := d
			cblas_dcopy(dim, vec_d, 1, vec_b, 1);

			printf("A := C! cond1(A) = %15.7e -> cond1(C) = %15.7e\n", cond_a, cond_c);
			
			// A has been changed!
			acculate_flag = 1;
		}
	
		old_norm_z = norm_z;

		// ||x||
		norm_x = cblas_dnrm2(dim, vec_x, 1);

		// r := vec_b - mat_a * x
		//cblas_dcopy(dim, vec_b, 1, vec_r, 1);
		//cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, -1.0, mat_a, dim, vec_x, 1, 1.0, vec_r, 1);
		accurate_dgemv(vec_r, mat_a, dim, dim, vec_x);
		cblas_daxpy(dim, -1.0, vec_b, 1, vec_r, 1);
		cblas_dscal(dim, -1.0, vec_r, 1);

		// z := r
		cblas_dcopy(dim, vec_r, 1, vec_z, 1);
	}

	// ret := x
	cblas_dcopy(dim, vec_x, 1, ret, 1);

	// �s�{�b�g�E�x�N�g������
	free(pivot);
	free(mat_lu);
	free(vec_z);
	free(vec_d);
	free(vec_r);
	free(vec_x);

	return itimes;
}

int main()
{
	lapack_int i, j, dim, itimes;
	lapack_int inc_vec_x, inc_vec_b;
	lapack_int *pivot, info;

	double *mat_a, *vec_b, *vec_x, *vec_ret;
	double *org_mat_a, *org_vec_b;
	double alpha, beta;
	double running_time;

	// ����������
	printf("Dim = "); scanf("%d", &dim);

	if(dim <= 0)
	{
		printf("Illigal dimenstion! (dim = %d)\n", dim);
		return EXIT_FAILURE;
	}

	// �ϐ�������
	mat_a = (double *)calloc(dim * dim, sizeof(double));
	org_mat_a = (double *)calloc(dim * dim, sizeof(double));
	vec_x = (double *)calloc(dim, sizeof(double));
	vec_b = (double *)calloc(dim, sizeof(double));
	org_vec_b = (double *)calloc(dim, sizeof(double));
	vec_ret = (double *)calloc(dim, sizeof(double));
	pivot = (lapack_int *)calloc(dim, sizeof(lapack_int));

	// mat_a��vec_x�ɒl����
	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
		{
			mat_a[i * dim + j] = (double)rand();
			//mat_a[i * dim + j] = (double)rand() / (double)RAND_MAX;
			//mat_a[i * dim + j] = 1.0 / (double)(i + j + 1);
			//if((i + j + 1) % 2 != 0)
			//	mat_a[i * dim + j] = -mat_a[i * dim + j];

			//mat_a[i * dim + j] = 1.0 / (double)(i + j + 1);
		}
		//mat_a[i * dim + i] += 2.0;
		//vec_x[i] = 1.0 / (double)(i + 1);
		vec_x[i] = (double)(i + 1);
	}

	// size(vec_x) == size(vec_b)
	inc_vec_x = inc_vec_b = 1;

	// vec_b := 1.0 * mat_a * vec_x + 0.0 * vec_b
	alpha = 1.0;
	beta = 0.0;
	cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, alpha, mat_a, dim, vec_x, inc_vec_x, beta, vec_b, inc_vec_b);

	// org_mat_a := mat_a
	// org_vec_b := vec_b
	cblas_dcopy(dim * dim, mat_a, 1, org_mat_a, 1);
	cblas_dcopy(dim, vec_b, 1, org_vec_b, 1);

	// �o��
/*	for(i = 0; i < dim; i++)
	{
		printf("[");
		for(j = 0; j < dim; j++)
			printf("%10.3f ", mat_a[i * dim + j]);
		printf("]  %10.3f = %10.3f\n", vec_x[i], vec_b[i]);
	}
*/
	// ���ږ@
	printf("Direct Methods:%d x %d\n", dim, dim);
	cblas_dcopy(dim, vec_b, 1, vec_ret, 1);
	itimes = 0;
	LAPACKE_dgesv(LAPACK_ROW_MAJOR, dim, 1, mat_a, dim, pivot, vec_ret, 1);

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

	// �������ǖ@
	printf("Iterative_ref Methods:%d x %d\n", dim, dim);
	cblas_dcopy(dim * dim, org_mat_a, 1, mat_a, 1);
	cblas_dcopy(dim, org_vec_b, 1, vec_b, 1);

	itimes = d_iterative_refinement(vec_ret, mat_a, vec_b, dim, 10, 1.0e-15, 0.0);
	//itimes = d_iterative_refinement_amm(vec_ret, mat_a, vec_b, dim, 10, 1.0e-15, 0.0);
	//itimes = d_iterative_refinement_amm(vec_ret, mat_a, vec_b, dim, 10, 0.0, 0.0);

	// print
	printf("iteration times = %d\n", itimes);
	printf("calculated x = \n");
/*	for(i = 0; i < dim; i++)
	{
		//printf("%3d -> %3d: ", i, pivot[i]);
		printf("%25.17e ", vec_ret[i]);
		printf("\n");
	}
*/
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
	
	// �̈�J��
	free(mat_a);
	free(org_mat_a);
	free(vec_x);
	free(vec_b);
	free(vec_ret);
	free(pivot);

	return EXIT_SUCCESS;
}
