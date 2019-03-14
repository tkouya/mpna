/*************************************************/
/* matmul.cc : 実行列×実行列
/* [Intel] icpc matmul.cc
/* [GCC  ] g++ matmul.cc
/*************************************************/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib> // EXIT_SUCCESS & EXIT_FAUILURE の定義

// matmul_gflops, byte_double_sqmt
#include "matmul_block.h"

// 時間計測用: get_secv, get_real_secv
#include "get_secv.h"

using namespace std;

// 正方行列×正方行列
// Row major方式
void matmul_simple(double ret[], double mat_a[], double mat_b[], int dim)
{
	int i, j, k, ij_index;

	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
		{
			ij_index = i * dim + j;
			ret[ij_index] = 0.0;
			for(k = 0; k < dim; k++)
				ret[ij_index] += mat_a[i * dim + k] * mat_b[k * dim + j];
		}
	}
}

// メイン関数
int main(int argc, char *argv[])
{
	int i, j, min_dim, max_dim, dim, iter, max_iter = 10, num_threads;
	double *mat_a, *mat_b, *mat_c;
	double stime, etime;

	if(argc < 3)
	{
		cout << "Usage: " << argv[0] << " [min. dimension]  [max.dimension]"<< endl;
		return EXIT_SUCCESS;
	}

	min_dim = atoi(argv[1]);
	max_dim = atoi(argv[2]);

	if(min_dim <= 0)
	{
		cout << "Illegal dimension! (min_dim = " << min_dim << ")" << endl;
		return EXIT_FAILURE;
	}

	// メインループ
	cout << setw(5) << "  dim :     SECONDS GFLOPS Mat.KB ||C||_F" << endl;
	for(dim = min_dim; dim <= max_dim; dim += 16)
	{

		// 変数初期化
		//mat_a = new double[dim * dim];
		mat_a = (double *)calloc(dim * dim, sizeof(double));
		mat_b = (double *)calloc(dim * dim, sizeof(double));
		mat_c = (double *)calloc(dim * dim, sizeof(double));

		// mat_aとmat_bに値入力
		for(i = 0; i < dim; i++)
		{
			for(j = 0; j < dim; j++)
			{
				mat_a[i * dim + j] = sqrt(5.0) * (double)(i + j + 1);
				mat_b[i * dim + j] = sqrt(3.0) * (double)(dim - (i + j));
			}
		}

		max_iter = 3; // 行列積を最低3回実行

		do
		{
			stime = get_real_secv();
			for(iter = 0; iter < max_iter; iter++)
				matmul_simple(mat_c, mat_a, mat_b, dim);
			etime = get_real_secv(); etime -= stime;

			if(etime >= 1.0) break; // 1秒以上になるまで繰り返し

			max_iter *= 2;
		} while(0);

		etime /= (double)max_iter; // 平均計算時間を導出
 
		// 出力
		cout << setw(5) << dim << " : " << setw(10) << setprecision(5) << etime << " " << matmul_gflops(etime, dim) << " " << byte_double_sqmat(dim) / 1024 << " " << normf_dmatrix_array(mat_c, dim, dim) << endl;

		// 変数消去
		free(mat_a);
		free(mat_b);
		free(mat_c);

	} // メインループ終了

	return EXIT_SUCCESS;
}
