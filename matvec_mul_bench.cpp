/*************************************************/
/* matvec_mul.cpp : 実行列×実ベクトル
/* icpt matvec_mul.cpp 
/*************************************************/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib> // define EXIT_SUCCESS & EXIT_FAUILURE

#include "get_secv.h"

using namespace std;

// GFlops
inline double matvec_mul_gflops(double comp_sec, int dim)
{
	return 2.0 * (double)dim * (double)dim / comp_sec / 1024.0 / 1024.0 / 1024.0;
}

// GB of Double prec. Square matrix
inline int byte_double_sqmat(int dim)
{
	return sizeof(double) * dim * dim;
}

// GB of Double prec. vector
inline int byte_double_vec(int dim)
{
	return sizeof(double) * dim;
}

// 正方行列×ベクトル
void matvec_mul_simple(double ret[], double mat[], double vec[], int dim)
{
	int i, j;

	for(i = 0; i < dim; i++)
	{
		ret[i] = 0.0;
		for(j = 0; j < dim; j++)
			ret[i] += mat[i * dim + j] * vec[j];
	}
}

// 正方行列×ベクトル
void matvec_mul_block(double ret[], double mat[], double vec[], int dim, int block_dim)
{
	int i, j;

	if(dim <= block_dim)
	{
		matvec_mul_simple(ret, mat, vec, dim);
		return;
	}

	// Part 1
	for(i = 0; i < block_dim; i++)
	{
		ret[i] = 0.0;
		for(j = 0; j < block_dim; j++)
			ret[i] += mat[i * dim + j] * vec[j];
	}

	// Part 2
	for(i = 0; i < block_dim; i++)
	{
		for(j = block_dim; j < dim; j++)
			ret[i] += mat[i * dim + j] * vec[j];
	}

	// Part 3
	for(i = block_dim; i < dim; i++)
	{
		ret[i] = 0.0;
		for(j = 0; j < block_dim; j++)
			ret[i] += mat[i * dim + j] * vec[j];
	}

	// Part 4
	for(i = block_dim; i < dim; i++)
	{
		for(j = block_dim; j < dim; j++)
			ret[i] += mat[i * dim + j] * vec[j];
	}
}

// 正方行列×ベクトル
void matvec_mul_simple_col(double ret[], double mat[], double vec[], int dim)
{
	int i, j;

	for(i = 0; i < dim; i++)
		ret[i] = 0.0;

	for(j = 0; j < dim; j++)
	{
		for(i = 0; i < dim; i++)
			ret[i] += mat[i + j * dim] * vec[j];
	}
}

int main(int argc, char *argv[])
{
	int i, j, min_dim, max_dim, dim, iter, max_iter = 10;
	double *mat_a, *vec_b, *vec_x;
	double stime, etime;

	// 次元数入力
//	cout << "DIM = ";
//	cin >> dim;

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
	for(dim = min_dim; dim <= max_dim; dim += 128)
	{

		// 変数初期化
		//mat_a = new double[dim * dim];
		mat_a = (double *)calloc(dim * dim, sizeof(double));
		vec_x = new double[dim];
		vec_b = new double[dim];

		// mat_aとvec_xに値入力
		for(i = 0; i < dim; i++)
		{
			for(j = 0; j < dim; j++)
			{
				mat_a[i * dim + j] = (double)(i + j + 1);
				if((i + j + 1) % 2 != 0)
					mat_a[i * dim + j] *= -1.0;
			}
			vec_x[i] = 1.0 / (double)(i + 1);
		}

		// 行列×ベクトル
		max_iter = 3;
		do
		{
			stime = get_secv();
			for(iter = 0; iter < max_iter; iter++)
				matvec_mul_block(vec_b, mat_a, vec_x, dim, 8192);
				//matvec_mul_simple_col(vec_b, mat_a, vec_x, dim);
				//matvec_mul_simple(vec_b, mat_a, vec_x, dim);
			etime = get_secv(); etime -= stime;
			if(etime >= 1.0) break;
			max_iter *= 8;
		} while(0);

		etime /= (double)max_iter;
 
		// 出力
		//cout << "Dimension     : " << dim << " * " << dim << endl;
		//cout << "Comp.Time(sec): " << setprecision(3) << etime << endl;
		//cout << "Gflops        : " << setprecision(3) << matvec_mul_gflops(etime, dim) << endl;
		cout << setw(5) << dim << " : " << matvec_mul_gflops(etime, dim) << " " << byte_double_vec(dim) / 1024 << endl;

		/*
		for(i = 0; i < dim; i++)
		{
			cout << " [ ";
			for(j = 0; j < dim; j++)
				cout << scientific << setprecision(3) << setw(10) << mat_a[i * dim + j] << " ";
			cout << " ] [ " << scientific << setprecision(3) << setw(10) << vec_x[i] << " ]   [ " << vec_b[i] << endl;
		}
		*/

		// 変数消去
		//delete mat_a;
		free(mat_a);
		delete vec_x;
		delete vec_b;

	} // end of mainloop

	return EXIT_SUCCESS;
}
