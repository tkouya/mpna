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

int main(int argc, char *argv[])
{
	int i, j, dim;

	double *mat_a, *vec_b, *vec_x;

	// 次元数入力
//	cout << "DIM = ";
//	cin >> dim;

	if(argc <= 1)
	{
		cout << "Usage: " << argv[0] << " [dimension]"<< endl;
		return EXIT_SUCCESS;
	}

	dim = atoi(argv[1]);

	if(dim <= 0)
	{
		cout << "Illegal dimension! (dim = " << dim << ")" << endl;
		return EXIT_FAILURE;
	}

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
	double stime = get_secv();
	matvec_mul_simple(vec_b, mat_a, vec_x, dim);
	double etime = get_secv(); etime -= stime;

	// 出力
	cout << "Dimension     : " << dim << " * " << dim << endl;
	cout << "Comp.Time(sec): " << setprecision(3) << etime << endl;
	cout << "Gflops        : " << setprecision(3) << matvec_mul_gflops(etime, dim) << endl;
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

	return EXIT_SUCCESS;
}
