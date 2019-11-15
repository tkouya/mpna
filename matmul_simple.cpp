//******************************************************************************
// matmul_simple.cpp : Matrix multiplication based on simple triple loop way
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
#include <cmath>

// matmul_gflops, byte_double_sqmt, and normf_dmatrix_array
#include "matmul_block.h"

// Time estimation: get_secv, get_real_secv
#include "get_secv.h"

using namespace std;

// square matrix multiplication (row-major order)
void matmul_simple(double ret[], double mat_a[], double mat_b[], int dim)
{
	int i, j, k, ij_index;

	for(i = 0; i < dim; i++)
	{
		for(j = 0; j < dim; j++)
		{
			ij_index = i * dim + j; // row-major order
			ret[ij_index] = 0.0;
			for(k = 0; k < dim; k++)
				ret[ij_index] += mat_a[i * dim + k] * mat_b[k * dim + j];
		}
	}
}

// main function
int main(int argc, char *argv[])
{
	int i, j, min_dim, max_dim, dim, iter, max_iter = 10, num_threads;
	double *mat_a, *mat_b, *mat_c;
	double stime, etime;

	if(argc < 3)
	{
		cout << "Usage: " << argv[0] << " [min. dimension]  [max.dimension]"<< endl;
		return 1;
	}

	min_dim = atoi(argv[1]);
	max_dim = atoi(argv[2]);

	if(min_dim <= 0)
	{
		cout << "Illegal dimension! (min_dim = " << min_dim << ")" << endl;
		return 1;
	}

	// main loop
	cout << setw(5) << "  dim :     SECONDS GFLOPS Mat.KB ||C||_F" << endl;
	for(dim = min_dim; dim <= max_dim; dim += 16)
	{

		// initialize matrices
		//mat_a = new double[dim * dim];
		mat_a = (double *)calloc(dim * dim, sizeof(double));
		mat_b = (double *)calloc(dim * dim, sizeof(double));
		mat_c = (double *)calloc(dim * dim, sizeof(double));

		// set values to mat_a and mat_b
		for(i = 0; i < dim; i++)
		{
			for(j = 0; j < dim; j++)
			{
				mat_a[i * dim + j] = sqrt(5.0) * (double)(i + j + 1);
				mat_b[i * dim + j] = sqrt(3.0) * (double)(dim - (i + j));
			}
		}

		max_iter = 3; // the number of trial time is 3 at least

		do
		{
			stime = get_real_secv();
			for(iter = 0; iter < max_iter; iter++)
				matmul_simple(mat_c, mat_a, mat_b, dim);
			etime = get_real_secv(); etime -= stime;

			if(etime >= 1.0) break; // repeat until total time is over 1 second

			max_iter *= 2;
		} while(0);

		etime /= (double)max_iter; // average time
 
		// output
		cout << setw(5) << dim << " : " << setw(10) << setprecision(5) << etime << " " << matmul_gflops(etime, dim) << " " << byte_double_sqmat(dim) / 1024 << " " << normf_dmatrix_array(mat_c, dim, dim) << endl;

		// delete matrices
		free(mat_a);
		free(mat_b);
		free(mat_c);

	} // end of main loop

	return 0;
}
