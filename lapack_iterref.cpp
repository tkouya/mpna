/*************************************************/
/* linear_eq_dgetrf.c : 連立一次方程式求解
/* icc linear_eq.c -L/usr/local/lib -llapacke -llapack -lcblas -lrefblas -L/opt/intel/lib/intel64 -lifcore
/*************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lapacke/lapacke.h"
#include "cblas.h"

// Iterative ref.
// ret[dim] := A[dim * dim]^(-1) * b[dim]
int d_iterative_refinement(double ret[], double mat_a[], double vec_b[], int dim, int max_itimes, double rtol, double atol)
{
	lapack_int info, itimes, i;
	lapack_int *pivot;
	double *vec_z, *vec_r, *vec_x, *mat_lu;
	double norm_z, norm_x;

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
	for(i = 0; i < dim; i++)
		vec_x[i] = 0.0;

	norm_x = cblas_dnrm2(dim, vec_x, 1);
	cblas_dcopy(dim, vec_b, 1, vec_r, 1);
	cblas_dcopy(dim, vec_r, 1, vec_z, 1);

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

		// 収束判定
		// ||(x + z) - x|| = ||z|| <= rtol * ||x|| + atol ?
		norm_z = cblas_dnrm2(dim, vec_z, 1);
		if(norm_z <= rtol * norm_x + atol)
			break;

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


int main()
{
	lapack_int i, j, dim, itimes;
	lapack_int inc_vec_x, inc_vec_b;
	lapack_int *pivot, info;

	double *mat_a, *vec_b, *vec_x, *vec_ret;
	double alpha, beta;
	double running_time;

	// 次元数入力
	printf("Dim = "); scanf("%d", &dim);

	if(dim <= 0)
	{
		printf("Illigal dimenstion! (dim = %d)\n", dim);
		return EXIT_FAILURE;
	}

	// 変数初期化
	mat_a = (double *)calloc(dim * dim, sizeof(double));
	vec_x = (double *)calloc(dim, sizeof(double));
	vec_b = (double *)calloc(dim, sizeof(double));
	vec_ret = (double *)calloc(dim, sizeof(double));

	// mat_aとvec_xに値入力
	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
		{
			//mat_a[i * dim + j] = (double)rand() / (double)RAND_MAX;
			mat_a[i * dim + j] = 1.0 / (double)(i + j + 1);
			if((i + j + 1) % 2 != 0)
				mat_a[i * dim + j] = -mat_a[i * dim + j];
		}
		mat_a[i * dim + i] += 2.0;
		vec_x[i] = 1.0 / (double)(i + 1);
	}

	// size(vec_x) == size(vec_b)
	inc_vec_x = inc_vec_b = 1;

	// vec_b := 1.0 * mat_a * vec_x + 0.0 * vec_b
	alpha = 1.0;
	beta = 0.0;
	cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, alpha, mat_a, dim, vec_x, inc_vec_x, beta, vec_b, inc_vec_b);

	// 出力
	for(i = 0; i < dim; i++)
	{
		printf("[");
		for(j = 0; j < dim; j++)
			printf("%10.3f ", mat_a[i * dim + j]);
		printf("]  %10.3f = %10.3f\n", vec_x[i], vec_b[i]);
	}

	// 反復改良法
	itimes = d_iterative_refinement(vec_ret, mat_a, vec_b, dim, 10, 1.0e-15, 0.0);

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
	printf("x - calculated x = \n");
	for(i = 0; i < dim; i++)
	{
		printf("%3d: ", i);
		printf("%10.2e ", fabs((vec_x[i] - vec_ret[i]) / vec_x[i]));
		printf("\n");
	}
	
	// 領域開放
	free(mat_a);
	free(vec_x);
	free(vec_b);
	free(vec_ret);

	return EXIT_SUCCESS;
}
