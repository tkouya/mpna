#ifndef __MATMUL_BLOCK_H
#define __MATMUL_BLOCK_H

#if defined (__cplusplus)
extern "C" {
#endif

// DVector 

typedef struct{
	double *element;
	long int dim;
} dvector;

typedef dvector *DVector;

#define get_dvector_i(vec, index) (*((vec)->element + (index)))

//double get_dvector_i(DVector vec, long int index)
//{
//	return *(vec->element + index);
//}

#define set_dvector_i(vec, index, val) (*((vec)->element + (index)) = (val))

//void set_dvector_i(DVector vec, long int index, double val)
//{
//	*(vec->element + index) = val;
//}

/*************************************************/
/* Vector Calculations for DVector               */
/*
DVector init_dvector(long int dimension)
void free_dvector(DVector vec)
void add_dvector(DVector c, DVector a, DVector b)
void add2_dvector(DVector c, DVector a)
void sub_dvector(DVector c, DVector a, DVector b)
void sub2_dvector(DVector c, DVector a)
void cmul_dvector(DVector c, double val, DVector a)
void cmul2_dvector(DVector c, double val)
void add_cmul_dvector(DVector c, DVector a, double val, DVector b)
double ip_dvector(DVector a, DVector b)
double norm1_dvector(DVector a)
double norm2_dvector(DVector a)
double normi_dvector(DVector a)
void subst_dvector(DVector c, DVector a)
*/
/*************************************************/

DVector init_dvector(long int dimension)
{
	DVector ret = NULL;
	long int i;

	if(dimension <= 0)
	{
		fprintf(stderr, "ERROR: init_dvector\n");
		return ret;
	}

	ret = (DVector)malloc(sizeof(dvector));
	if(ret == NULL)
		return ret;

	ret->element = (double *)calloc(sizeof(double), dimension);
	if(ret->element == NULL)
		return ret;

	/* All 0 */
	for(i = 0; i < dimension; i++)
		*(ret->element + i) = 0.0;

	ret->dim = dimension;

	return ret;
}

void free_dvector(DVector vec)
{
	if(vec == NULL)
		return;

	if(vec->element != NULL)
		free(vec->element);

	free(vec);
}

/* c = a + b */
void add_dvector(DVector c, DVector a, DVector b)
{
	long int i;

	if((a->dim != b->dim) || (c->dim != a->dim) || (c->dim != b->dim))
	{
		fprintf(stderr, "ERROR: add_dvector\n");
		return;
	}

	for(i = 0; i < c->dim; i++)
		set_dvector_i(c, i, get_dvector_i(a, i) + get_dvector_i(b, i));

}

/* c += a */
void add2_dvector(DVector c, DVector a)
{
	long int i;

	if(c->dim != a->dim)
	{
		fprintf(stderr, "ERROR: add2_dvector\n");
		return;
	}

	for(i = 0; i < c->dim; i++)
		set_dvector_i(c, i, get_dvector_i(c, i) + get_dvector_i(a, i));

}

/* c = a - b */
void sub_dvector(DVector c, DVector a, DVector b)
{
	long int i;

	if((a->dim != b->dim) || (c->dim != a->dim) || (c->dim != b->dim))
	{
		fprintf(stderr, "ERROR: sub_dvector\n");
		return;
	}

	for(i = 0; i < c->dim; i++)
		set_dvector_i(c, i, get_dvector_i(a, i) - get_dvector_i(b, i));

}

/* c -= a */
void sub2_dvector(DVector c, DVector a)
{
	long int i;

	if(c->dim != a->dim)
	{
		fprintf(stderr, "ERROR: sub2_dvector\n");
		return;
	}

	for(i = 0; i < c->dim; i++)
		set_dvector_i(c, i, get_dvector_i(c, i) - get_dvector_i(a, i));

}

/* c = val * a */
void cmul_dvector(DVector c, double val, DVector a)
{
	long int i;

	if(c->dim != a->dim)
	{
		fprintf(stderr, "ERROR: cmul_dvector\n");
		return;
	}

	for(i = 0; i < c->dim; i++)
		set_dvector_i(c, i, val * get_dvector_i(a, i));

}

/* c *= val */
void cmul2_dvector(DVector c, double val)
{
	long int i;

	for(i = 0; i < c->dim; i++)
		set_dvector_i(c, i, val * get_dvector_i(c, i));

}

/* c = a + val * b */
void add_cmul_dvector(DVector c, DVector a, double val, DVector b)
{
	long int i;

	if((a->dim != b->dim) || (c->dim != a->dim) || (c->dim != b->dim))
	{
		fprintf(stderr, "ERROR: add_cmul_dvector\n");
		return;
	}

	for(i = 0; i < c->dim; i++)
		set_dvector_i(c, i, get_dvector_i(a, i) + val * get_dvector_i(b, i));

}

/* (a, b) */
double ip_dvector(DVector a, DVector b)
{
	double tmp = 0.0;
	long int i;

	if(a->dim != b->dim)
	{
		fprintf(stderr, "ERROR: ip_dvector\n");
		return 0;
	}

	for(i = 0; i < a->dim; i++)
		tmp += get_dvector_i(a, i) * get_dvector_i(b, i);

	return tmp;
}

/* ||a||_1 */
double norm1_dvector(DVector a)
{
	double ret = 0.0;
	long int i;

	for(i = 0; i < a->dim; i++)
		ret += fabs(get_dvector_i(a, i));

	return ret;
}

/* ||a||_2 */
double norm2_dvector(DVector a)
{
	double ret = 0.0;
	long int i;

	for(i = 0; i < a->dim; i++)
		ret += get_dvector_i(a, i) * get_dvector_i(a, i);

	return sqrt(ret);
}

/* ||a||_infty */
double normi_dvector(DVector a)
{
	double ret, tmp;
	long int i;

	ret = fabs(get_dvector_i(a, 0));
	for(i = 1; i < a->dim; i++)
	{
		tmp = fabs(get_dvector_i(a, i));
		if(ret < tmp)
			ret = tmp;
	}

	return ret;
}

/* c := a */
void subst_dvector(DVector c, DVector a)
{
	long int i;

	for(i = 0; i < a->dim; i++)
		set_dvector_i(c, i, get_dvector_i(a, i));
}

/* c := 0 */
void set0_dvector(DVector c)
{
	long int i;

	for(i = 0; i < c->dim; i++)
		set_dvector_i(c, i, (double)0);
}

/* append 2005.07/12 */
/*
	ret(index_start) = src(src_index_start)
	 ...
	ret(index_end  ) = src(src_index_end)
*/
void copy_dvector_ij(DVector ret, long int index_start, long int index_end, DVector src, long int src_index_start, long int src_index_end)
{
	long int i, itmp;

	if((src_index_end - src_index_start) != (index_end - index_start))
	{
		fprintf(stderr, "Invalid index!(copy_dvector_ij)\n");
		return;
	}

	for(i = 0; i <= (index_end - index_start); i++)
	{
		set_dvector_i(ret, index_start + i, get_dvector_i(src, src_index_start + i));
//		printf("%d <----------------------------------> %d\n", index_start + i, src_index_start + i);
	}
}

// DMatrix

typedef struct{
	double *element;
	long int row_dim, col_dim;
} dmatrix;

typedef dmatrix *DMatrix;

#define get_dmatrix_ij(mat, i, j) ( *((mat)->element + (i) * (mat)->col_dim + (j)) )
#define set_dmatrix_ij(mat, i, j, val) ( *((mat)->element + (i) * (mat)->col_dim + (j)) = (val) )


/*************************************************/
/* Matrix Caluculations for DMatrix              */
/*
DMatrix init_dmatrix(long int row_dimension, long int col_dimension)
void free_dmatrix(DMatrix mat)
double normf_dmatrix(DMatrix mat)
double normi_dmatrix(DMatrix mat)
double norm1_dmatrix(DMatrix mat)
void add_dmatrix(DMatrix c, DMatrix a, DMatrix b);
void sub_dmatrix(DMatrix c, DMatrix a, DMatrix b);
void mul_dmatrix(DMatrix c, DMatrix a, DMatrix b);
void transpose_dmatrix(DMatrix c, DMatrix a);
void mul_dmatrix_dvec(DVector v, DMatrix a, DVector vb)
void mul_dmatrixt_dvec(DVector v, DMatrix a, DVector vb)
void inv_dmatrix(DMatrix a);
void subst_dmatrux(DMatrix c, DMatrix a);
*/
/*************************************************/
DMatrix init_dmatrix(long int row_dimension, long int col_dimension)
{
	DMatrix ret = NULL;
	long int i, j;

	if(row_dimension <= 0 || col_dimension <= 0)
	{
		fprintf(stderr, "ERROR: init_dmatrix\n");
		return ret;
	}

	ret = (DMatrix)malloc(sizeof(dmatrix));
	if(ret == NULL)
		return ret;

	ret->element = (double *)calloc(sizeof(double), row_dimension * col_dimension);
	if(ret->element == NULL)
		return ret;

	/* All 0 */
	for(i = 0; i < row_dimension; i++)
		for(j = 0; j < col_dimension; j++)
			*(ret->element + i * col_dimension + j) = (double)0.0;

	ret->row_dim = row_dimension;
	ret->col_dim = col_dimension;

	return ret;
}

void free_dmatrix(DMatrix mat)
{
	if(mat == NULL)
		return;

	if(mat->element != NULL)
		free(mat->element);

	free(mat);
}

/*************************************************/
/* Matrix Caluculations for DMatrix              */
/*
double normf_dmatrix(DMatrix mat)
double normi_dmatrix(DMatrix mat)
double norm1_dmatrix(DMatrix mat)
void add_dmatrix(DMatrix c, DMatrix a, DMatrix b);
void sub_dmatrix(DMatrix c, DMatrix a, DMatrix b);
void mul_dmatrix(DMatrix c, DMatrix a, DMatrix b);
void transpose_dmatrix(DMatrix c, DMatrix a);
void mul_dmatrix_dvec(DVector v, DMatrix a, DVector vb)
void mul_dmatrixt_dvec(DVector v, DMatrix a, DVector vb)
void inv_dmatrix(DMatrix a);
void subst_dmatrux(DMatrix c, DMatrix a);
*/
/*************************************************/
/* Frobenius Norm of Matrix */
double normf_dmatrix(DMatrix mat)
{
	long int i, j;
	double ret;

	ret = 0.0;
	for(i = 0; i < mat->row_dim; i++)
		for(j = 0; j < mat->col_dim; j++)
			ret += (get_dmatrix_ij(mat, i, j) * get_dmatrix_ij(mat, i, j));

	ret = sqrt(ret);

	return ret;
}

/* Frobenius Norm of Matrix: array type */
double normf_dmatrix_array(double mat[], int row_dim, int col_dim)
{
	int i, j;
	double ret;

	ret = 0.0;
	for(i = 0; i < row_dim; i++)
		for(j = 0; j < col_dim; j++)
			ret += mat[i * col_dim + j] * mat[i * col_dim + j];

	ret = sqrt(ret);

	return ret;
}

/* Infinity Norm of Matrix */
double normi_dmatrix(DMatrix mat)
{
	long int i, j;
	double ret, sum;

	ret = 0.0;
	for(i = 0; i < mat->row_dim; i++)
	{
		sum = 0.0;
		for(j = 0; j < mat->col_dim; j++)
			sum += fabs(get_dmatrix_ij(mat, i, j));
		if(ret < sum)
			ret = sum;
	}

	return ret;
}

/* 1 Norm of Matrix */
double norm1_dmatrix(DMatrix mat)
{
	long int i, j;
	double ret, sum;

	ret = 0.0;
	for(j = 0; j < mat->col_dim; j++)
	{
		sum = 0.0;
		for(i = 0; i < mat->row_dim; i++)
			sum += fabs(get_dmatrix_ij(mat, i, j));
		if(ret < sum)
			ret = sum;
	}

	return ret;
}

/* c = a + b */
void add_dmatrix(DMatrix c, DMatrix a, DMatrix b)
{
	long int i, j, row_dim, col_dim;

	/* check row_dim */
	if((a->row_dim != b->row_dim) || (b->row_dim != c->row_dim) || (c->row_dim != a->row_dim))
	{
		fprintf(stderr, "ERROR: add_dmatrix\n");
		return;
	}
	row_dim = c->row_dim;

	/* check col_dim */
	if((a->col_dim != b->col_dim) || (b->col_dim != c->col_dim) || (c->col_dim != a->col_dim))
	{
		fprintf(stderr, "ERROR: add_dmatrix\n");
		return;
	}
	col_dim = c->col_dim;

	for(i = 0; i < row_dim; i++)
	{
		for(j = 0; j < col_dim; j++)
			set_dmatrix_ij(c, i, j, get_dmatrix_ij(a, i, j) + get_dmatrix_ij(b, i, j));
	}
}

/* c = a - b */
void sub_dmatrix(DMatrix c, DMatrix a, DMatrix b)
{
	long int i, j, row_dim, col_dim;

	/* check row_dim */
	if((a->row_dim != b->row_dim) || (b->row_dim != c->row_dim) || (c->row_dim != a->row_dim))
	{
		fprintf(stderr, "ERROR: sub_dmatrix\n");
		return;
	}
	row_dim = c->row_dim;

	/* check col_dim */
	if((a->col_dim != b->col_dim) || (b->col_dim != c->col_dim) || (c->col_dim != a->col_dim))
	{
		fprintf(stderr, "ERROR: sub_dmatrix\n");
		return;
	}
	col_dim = c->col_dim;

	for(i = 0; i < row_dim; i++)
	{
		for(j = 0; j < col_dim; j++)
			set_dmatrix_ij(c, i, j, get_dmatrix_ij(a, i, j) - get_dmatrix_ij(b, i, j));
	}
}

/* c = sc * a */
void cmul_dmatrix(DMatrix c, double sc, DMatrix a)
{
	long int i, j, row_dim, col_dim;

	/* check row_dim */
	if(a->row_dim != c->row_dim)
	{
		fprintf(stderr, "ERROR: cmul_dmatrix\n");
		return;
	}
	row_dim = c->row_dim;

	/* check col_dim */
	if(a->col_dim != c->col_dim)
	{
		fprintf(stderr, "ERROR: cmul_dmatrix\n");
		return;
	}
	col_dim = c->col_dim;

	for(i = 0; i < row_dim; i++)
	{
		for(j = 0; j < col_dim; j++)
			set_dmatrix_ij(c, i, j, sc * get_dmatrix_ij(a, i, j));
	}
}

/* c = a * b */
void mul_dmatrix(DMatrix c, DMatrix a, DMatrix b)
{
	long int i, j, k;
	double tmp;

	/* dimension check */
	if((c->row_dim != a->row_dim) || (c->col_dim != b->row_dim) || (a->col_dim != b->row_dim))
	{
		fprintf(stderr, "ERROR: mul_dmatrix\n");
		return;
	}

	for(i = 0; i < c->row_dim; i++)
	{
		for(j = 0; j < c->col_dim; j++)
		{
			tmp = 0.0;
			for(k = 0; k < a->col_dim; k++)
				tmp += get_dmatrix_ij(a, i, k) * get_dmatrix_ij(b, k, j);
			set_dmatrix_ij(c, i, j, tmp);
		}
	}
}

/* c = a^T */
void transpose_dmatrix(DMatrix c, DMatrix a)
{
	long int i, j;

	/* Check Dimentions */
	if((c->row_dim != a->col_dim) || (c->col_dim != a->row_dim))
	{
		fprintf(stderr, "ERROR: transpose_dmatrix\n");
		return;
	}
	
	for(i = 0; i < c->row_dim; i++)
	{
		for(j = 0; j < c->col_dim; j++)
			set_dmatrix_ij(c, i, j, get_dmatrix_ij(a, j, i));
	}
}

/* c := a */
void subst_dmatrix(DMatrix c, DMatrix a)
{
	long int i, j;

	if((c->row_dim != a->row_dim) || (c->col_dim != a->col_dim))
	{
		fprintf(stderr, "ERROR: subst_dmatrix\n");
		return;
	}

	for(i = 0; i < a->row_dim; i++)
	{
		for(j = 0; j < a->col_dim; j++)
		{
			set_dmatrix_ij(c, i, j, get_dmatrix_ij(a, i, j));
		}
	}
}

/* c := 0 */
void set0_dmatrix(DMatrix c)
{
	long int i, j;

	for(i = 0; i < c->row_dim; i++)
	{
		for(j = 0; j < c->col_dim; j++)
			set_dmatrix_ij(c, i, j, 0.0);
	}
}

/* c := I */
void setI_dmatrix(DMatrix c)
{
	long int i, j;

	for(i = 0; i < c->row_dim; i++)
	{
		for(j = 0; j < c->col_dim; j++)
			set_dmatrix_ij(c, i, j, 0.0);
		if(i < c->col_dim)
			set_dmatrix_ij(c, i, i, 1.0);
	}
}

/* v = a * vb */
void mul_dmatrix_dvec(DVector v, DMatrix a, DVector vb)
{
	long int i, j;
	double tmp;

	/* Check Dimension */
	if((v->dim != a->row_dim) || (vb->dim != a->col_dim))
	{
		fprintf(stderr, "ERROR: mul_dmatrix_dvec\n");
		return;
	}
	
	for(i = 0; i < a->row_dim; i++)
	{
		tmp = 0.0;
		for(j = 0; j < a->col_dim; j++)
			tmp += get_dmatrix_ij(a, i, j) * get_dvector_i(vb, j);
		set_dvector_i(v, i, tmp);
	}
}

/* v = a^T * vb */
void mul_dmatrixt_dvec(DVector v, DMatrix a, DVector vb)
{
	long int i, j;
	double tmp;

	/* Check Dimension */
	if((v->dim != a->col_dim) || (vb->dim != a->row_dim))
	{
		fprintf(stderr, "ERROR: mul_dmatrixt_dvec\n");
		return;
	}
	
	for(i = 0; i < a->col_dim; i++)
	{
		tmp = 0.0;
		for(j = 0; j < a->row_dim; j++)
			tmp += get_dmatrix_ij(a, j, i) * get_dvector_i(vb, j);
		set_dvector_i(v, i, tmp);
	}
}

// partial add
// *_index[0] = start_row_index
// *_index[1] = end_row_index
// *_index[2] = start_col_index
// *_index[3] = end_col_index
void add_dmatrix_partial(DMatrix ret, long int ret_index[4], DMatrix mat_a, long int mat_a_index[4], DMatrix mat_b, long int mat_b_index[4])
{
	long int i, j, ret_i, ret_j, a_i, a_j, b_i, b_j;
	long int imax, jmax;
	double tmp_val;

	imax = ret_index[1] - ret_index[0];
	jmax = ret_index[3] - ret_index[2];

	for(i = 0; i < imax; i++)
	{
		ret_i = ret_index[0] + i;
		a_i = mat_a_index[0] + i;
		b_i = mat_b_index[0] + i;
		//printf("i: %ld %ld %ld\n", ret_i, a_i, b_i);
		for(j = 0; j < jmax; j++)
		{
			ret_j = ret_index[2] + j;
			a_j = mat_a_index[2] + j;
			b_j = mat_b_index[2] + j;
			//printf("j: %ld %ld %ld\n", ret_j, a_j, b_j);

			tmp_val = get_dmatrix_ij(mat_a, a_i, a_j) + get_dmatrix_ij(mat_b, b_i, b_j);
			set_dmatrix_ij(ret, ret_i, ret_j, tmp_val);
		}
	}
}

// partial sub
// *_index[0] = start_row_index
// *_index[1] = end_row_index
// *_index[2] = start_col_index
// *_index[3] = end_col_index
void sub_dmatrix_partial(DMatrix ret, long int ret_index[4], DMatrix mat_a, long int mat_a_index[4], DMatrix mat_b, long int mat_b_index[4])
{
	long int i, j, ret_i, ret_j, a_i, a_j, b_i, b_j;
	long int imax, jmax;
	double tmp_val;

	imax = ret_index[1] - ret_index[0];
	jmax = ret_index[3] - ret_index[2];

	for(i = 0; i < imax; i++)
	{
		ret_i = ret_index[0] + i;
		a_i = mat_a_index[0] + i;
		b_i = mat_b_index[0] + i;
		for(j = 0; j < jmax; j++)
		{
			ret_j = ret_index[2] + j;
			a_j = mat_a_index[2] + j;
			b_j = mat_b_index[2] + j;

			tmp_val = get_dmatrix_ij(mat_a, a_i, a_j) - get_dmatrix_ij(mat_b, b_i, b_j);
			set_dmatrix_ij(ret, ret_i, ret_j, tmp_val);
		}
	}
}

// partial set
// *_index[0] = start_row_index
// *_index[1] = end_row_index
// *_index[2] = start_col_index
// *_index[3] = end_col_index
void subst_dmatrix_partial(DMatrix ret, long int ret_index[4], DMatrix mat_a, long int mat_a_index[4])
{
	long int i, j, ret_i, ret_j, a_i, a_j;
	long int imax, jmax;

	imax = ret_index[1] - ret_index[0];
	jmax = ret_index[3] - ret_index[2];

	for(i = 0; i < imax; i++)
	{
		ret_i = ret_index[0] + i;
		a_i = mat_a_index[0] + i;
		for(j = 0; j < jmax; j++)
		{
			ret_j = ret_index[2] + j;
			a_j = mat_a_index[2] + j;
			set_dmatrix_ij(ret, ret_i, ret_j, get_dmatrix_ij(mat_a, a_i, a_j));
		}
	}
}

// partial set
// *_index[0] = start_row_index
// *_index[1] = end_row_index
// *_index[2] = start_col_index
// *_index[3] = end_col_index
void neg_dmatrix_partial(DMatrix ret, long int ret_index[4], DMatrix mat_a, long int mat_a_index[4])
{
	long int i, j, ret_i, ret_j, a_i, a_j;
	long int imax, jmax;

	imax = ret_index[1] - ret_index[0];
	jmax = ret_index[3] - ret_index[2];

	for(i = 0; i < imax; i++)
	{
		ret_i = ret_index[0] + i;
		a_i = mat_a_index[0] + i;
		for(j = 0; j < jmax; j++)
		{
			ret_j = ret_index[2] + j;
			a_j = mat_a_index[2] + j;
			set_dmatrix_ij(ret, ret_i, ret_j, -get_dmatrix_ij(mat_a, a_i, a_j));
		}
	}
}

// partial set
// *_index[0] = start_row_index
// *_index[1] = end_row_index
// *_index[2] = start_col_index
// *_index[3] = end_col_index
void subst_dmatrix_partial_checked(DMatrix ret, long int ret_index[4], DMatrix mat_a, long int mat_a_index[4])
{
	long int i, j, ret_i, ret_j, a_i, a_j;
	long int imax, jmax;

	imax = ret_index[1] - ret_index[0];
	jmax = ret_index[3] - ret_index[2];

	for(i = 0; i < imax; i++)
	{
		ret_i = ret_index[0] + i;
		a_i = mat_a_index[0] + i;
		if((ret_i >= 0) && (ret_i < ret->row_dim))
		{
			for(j = 0; j < jmax; j++)
			{
				ret_j = ret_index[2] + j;
				a_j = mat_a_index[2] + j;
				if((ret_j >= 0) && (ret_j < ret->col_dim))
				{
					if((a_i >= 0) && (a_i < mat_a->row_dim) && (a_j >= 0) && (a_j < mat_a->col_dim))
						set_dmatrix_ij(ret, ret_i, ret_j, get_dmatrix_ij(mat_a, a_i, a_j));
					else
						set_dmatrix_ij(ret, ret_i, ret_j, 0.0); // Padding
					//printf("Warning: ret_index = %d, %d, %d, %d\n", ret_i, ret_j, a_i, a_j);
				}
			}
		}
	}
}

// log2(x)
double mylog2(double x);

// Padding to 2-powered dimensional matrix
DMatrix init_static_padding_dmatrix_strassen(DMatrix orig_mat);

// Dynamic Padding to even dimensional matrix
DMatrix init_dynamic_padding_dmatrix_strassen(DMatrix orig_mat);

// Strassen's Algorithm with static padding
void mul_dmatrix_strassen_odd_padding(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// Fit dimension to be multiple of min_dim
void mul_dmatrix_strassen(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// Strassen's Algorithm
void mul_dmatrix_strassen_odd_peeling(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// Strassen's Algorithm
void mul_dmatrix_strassen_even(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// Padding to 2-powered dimensional matrix
DMatrix init_static_padding_dmatrix_winograd(DMatrix orig_mat);

// Dynamic Padding to even dimensional matrix
DMatrix init_dynamic_padding_dmatrix_winograd(DMatrix orig_mat);

// Winograd's Algorithm with static padding
void mul_dmatrix_winograd_odd_padding(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// Fit dimension to be multiple of min_dim
void mul_dmatrix_winograd(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// Winograd's Algorithm
void mul_dmatrix_winograd_odd_peeling(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// Winograd's Algorithm
void mul_dmatrix_winograd_even(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// Block matrix multiplicaiton
void mul_dmatrix_blocked(DMatrix ret, DMatrix mat_a, DMatrix mat_b, long int min_dim);

// log2(x)
inline double mylog2(double x)
{
	return log10(x) / log10(2.0);
}

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

#if defined (__cplusplus)
}
#endif

#endif // __MATMUL_BLOCK_H
