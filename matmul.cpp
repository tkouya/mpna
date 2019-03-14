/*************************************************/
/* matmul.cc : 実行列×実行列
/* [Intel] icpc matmul.cc
/* [GCC  ] g++ matmul.cc
/*************************************************/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib> // define EXIT_SUCCESS & EXIT_FAUILURE

#include "get_secv.h"

using namespace std;

// GFlops
inline double matmul_gflops(double comp_sec, int dim)
{
	return 2.0 * (double)dim * (double)dim * (double)dim / comp_sec / 1024.0 / 1024.0 / 1024.0;
}

// GB of Double prec. Square matrix
inline int byte_double_sqmat(int dim)
{
	return sizeof(double) * dim * dim;
}

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

int main(int argc, char *argv[])
{
	int i, j, min_dim, max_dim, dim, iter, max_iter = 10, num_threads;
	double *mat_a, *mat_b, *mat_c;
	double stime, etime;

	// 次元数入力
//	cout << "DIM = ";
//	cin >> dim;

#ifdef _OPENMP
	cout << "num_threads: ";
	cin >> num_threads;

	omp_set_num_threads(num_threads);
#endif // _OPENMP

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

	// mainloop
	cout << setw(5) << "  dim :     SECONDS GFLOPS Mat.KB" << endl;
	//for(dim = min_dim; dim <= max_dim; dim += 128)
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
				mat_a[i * dim + j] = (double)(i + j + 1);
				if((i + j + 1) % 2 != 0)
					mat_a[i * dim + j] *= -1.0;

				mat_b[i * dim + j] = (double)(i + j + 1);
				if((i + j + 1) % 2 != 0)
					mat_b[i * dim + j] *= -1.0;
			}
		}

		// 行列×行列
		max_iter = 3;
		do
		{
			stime = get_secv();
			for(iter = 0; iter < max_iter; iter++)
				matmul_simple(mat_c, mat_a, mat_b, dim);
			etime = get_secv(); etime -= stime;

			if(etime >= 1.0) break;
			max_iter *= 2;
		} while(0);

		etime /= (double)max_iter;
 
		// 出力
		//cout << "Dimension     : " << dim << " * " << dim << endl;
		//cout << "Comp.Time(sec): " << setprecision(3) << etime << endl;
		//cout << "Gflops        : " << setprecision(3) << matmul_gflops(etime, dim) << endl;
		cout << setw(5) << dim << " : " << setw(10) << setprecision(5) << etime << " " << matmul_gflops(etime, dim) << " " << byte_double_sqmat(dim) / 1024 << endl;

		/*
		for(i = 0; i < dim; i++)
		{
			cout << " [ ";
			for(j = 0; j < dim; j++)
				cout << scientific << setprecision(3) << setw(10) << mat_a[i * dim + j] << " ";
			cout << "] " << endl;
		}
		*/

		// 変数消去
		//delete mat_a;
		free(mat_a);
		free(mat_b);
		free(mat_c);

	} // end of mainloop

	return EXIT_SUCCESS;
}
