#include <iostream>
#include <cmath>
#include <cstdio>
#include <qd/inline.h>
#include <qd/qd_real.h>

#include "mpfr.h"

using namespace std;

//using namespace qd;

// From QD package
#define _BOD_DOUBLE_SPLITTER		(134217729.0) // = 2^27 + 1
#define _BOD_DOUBLE_SPLIT_THRESH	(6.696928794914170755927656556625e+299) // = 2^996

// Original
//#define _BOD_FLOAT_SPLITTER		((float)(8193.0)) // = 2^13 + 1
#define _BOD_FLOAT_SPLITTER		((float)(4097.0)) // = 2^12 + 1
#define _BOD_FLOAT_SPLIT_THRESH	((float)(1.6615349947311448411297588253504e+35)) // = 2^117


// double to float[2] = {hi, lo}
inline void split_d2d(double &ret_high, double &ret_low, double val)
{
	double tmp;
	
	if((val > _BOD_DOUBLE_SPLIT_THRESH) || (val < _BOD_DOUBLE_SPLIT_THRESH))
	{
		val *= 3.7252902984619140625E-09; // * 2^(-28)
		tmp = _BOD_DOUBLE_SPLITTER * val;
		ret_high = tmp - (tmp - val);
		ret_low  = val - ret_high;
		ret_high *= 268435456.0; // * 2^28
		ret_low  *= 268435456.0; // * 2^28
	}
	else
	{
		tmp = _BOD_DOUBLE_SPLITTER * val;
		ret_high = (float)(tmp - (tmp - val));
		ret_low  = (float)(val - ret_high);
	}
}

// float[2] = {hi, lo} to double
inline double merge_f2d(float val_high, float val_low)
{
	return (double)((double)val_high + (double)val_low);
}

// float to float[2] = {hi, lo}
inline void split_f2f(float &ret_high, float &ret_low, float val)
{
	float tmp;
	
	if((val > _BOD_FLOAT_SPLIT_THRESH) || (val < _BOD_FLOAT_SPLIT_THRESH))
	{
//		val *= 6.103515625E-5; // 0.00006103515625 = 2^(-14)
		val *= 1.220703125E-4; // 0.0001220703125 = 2^(-13)
		tmp = _BOD_FLOAT_SPLITTER * val;
		ret_high  = (float)(tmp - (tmp - val));
		ret_low   = (float)(val - ret_high);
		//*ret_high  *= 16384.0; // * 2^14
		//*ret_low   *= 16384.0; // * 2^14
		ret_high  *= 8192.0; // * 2^13
		ret_low   *= 8192.0; // * 2^13
	}
	else
	{
		tmp = _BOD_FLOAT_SPLITTER * val;
		ret_high = (float)(tmp - (tmp - val));
		ret_low  = (float)(val - ret_high);
	}
}

// ----------------------------
// Double precision
// ----------------------------

// Quick-Two-Sum : assume |a| > |b|
// s := quick_two_sum(a, b, &error)
inline float fquick_two_sum(float a, float b, float &error)
{
	float s;

/*	if(b == 0.0)
	{
		*error = (float)0.0;
		return (a + b);
	}
*/
	s = a + b;
	error = b - (s - a);

	return s;
}

// float[2] = {hi, lo} to double
inline double merge_f2d_normalize(float val_high, float val_low)
{
	double ret;
	float error = 0.0;

	ret = (double)fquick_two_sum(val_high, val_low, error);

	return (double)(ret + (double)error);
}

inline double dquick_two_sum(double a, double b, double &error)
{
	double s;

	s = a + b;
	error = b - (s - a);

	return s;
}

#include "dd_print.h"

// inner product
double ddotp(double a[], double b[], int dim)
{
	double ret = 0;

	for(int i = 0; i < dim; i++)
		ret = fma(a[i], b[i], ret);

//		ret += a[i] * b[i];

}

// dotxblas
double ddotp_xblas(double a[], double b[], int dim)
{
	double s, t, h, r, s1, s2, t1, t2;

	s = 0; t = 0;

	for(int i = 0; i < dim; i++)
	{
		h = qd::two_prod(a[i], b[i], r);
		s1 = qd::two_sum(s, h, s2);
		t1 = qd::two_sum(t, r, t2);
		s2 += t1;
		t1 = qd::quick_two_sum(s1, s2, s2);
		t2 += s2;
		s = qd::quick_two_sum(t1, t2, t);
	}

	return s;
}

// mydotp : x^T * y
void mpfr_mydotp_d(mpfr_t ret, double x[], double y[], int dim)
{
	int i;
	mpfr_t x_element, y_element;

	mpfr_init2(x_element, mpfr_get_prec(ret));
	mpfr_init2(y_element, mpfr_get_prec(ret));

	mpfr_set_ui(ret, 0UL, MPFR_RNDN);

	for(i = 0; i < dim; i++)
	{
		mpfr_set_d(x_element, x[i], MPFR_RNDN);
		mpfr_set_d(y_element, y[i], MPFR_RNDN);

		mpfr_fma(ret, x_element, y_element, ret, MPFR_RNDN);
	}

	mpfr_clear(x_element);
	mpfr_clear(y_element);

	return;
}

// mydotp : x^T * y
double d_mpfrdotp_d(double x[], double y[], int dim)
{
	double dret;
	mpfr_t ret;

	mpfr_init2(ret, 512);

	mpfr_mydotp_d(ret, x, y, dim);

	dret = mpfr_get_d(ret, MPFR_RNDN);

	mpfr_clear(ret);

	return dret;
}

// mydotp : x^T * y
double d_mpfr_sum_d(double x[], int dim)
{
	double dret;
	mpfr_t ret;

	mpfr_init2(ret, 512);

	mpfr_set_d(ret, x[0], MPFR_RNDN);
	for(int i = 1; i < dim; i++)
		mpfr_add_d(ret, ret, x[i], MPFR_RNDN);

	dret = mpfr_get_d(ret, MPFR_RNDN);
	mpfr_clear(ret);

	return dret;
}

// mydotp : x^T * y
double relerr_ddotp(double x[], double y[], int dim)
{
	int i;
	double ddotp_val, relerr;
	mpfr_t x_element, y_element, ret, tmp;

	ddotp_val = ddotp(x, y, dim);

	mpfr_init2(ret, 512);

	mpfr_init2(tmp, mpfr_get_prec(ret));
	mpfr_init2(x_element, mpfr_get_prec(ret));
	mpfr_init2(y_element, mpfr_get_prec(ret));

	mpfr_set_ui(ret, 0UL, MPFR_RNDN);

	for(i = 0; i < dim; i++)
	{
		mpfr_set_d(x_element, x[i], MPFR_RNDN);
		mpfr_set_d(y_element, y[i], MPFR_RNDN);

		mpfr_fma(ret, x_element, y_element, ret, MPFR_RNDN);
	}

	mpfr_set_d(tmp, ddotp_val, MPFR_RNDN);
	mpfr_reldiff(ret, ret, tmp, MPFR_RNDN);
	relerr = mpfr_get_d(ret, MPFR_RNDN);

	mpfr_clear(tmp);
	mpfr_clear(x_element);
	mpfr_clear(y_element);
	mpfr_clear(ret);

	return relerr;
}

// error-free splitting
double split_dvec(double high_vec[], double low_vec[], const double org_vec[], int dim)
{
	int i;
	mpfr_t exact_val;
	double abs_a, beta, v, sigma;

	beta = ceil((log2((double)(dim + 1)) + 53.0) / 2.0);

	v = fabs(org_vec[0]);
	for(i = 1; i < dim; i++)
	{
		abs_a = fabs(org_vec[i]);
		if(abs_a > v)
			v = abs_a;
	}
	v = ceil(log2(v));

	sigma = pow(2.0, beta + v);
//	printf("sigma, tau = %10.3e, %10.3e\n", sigma, tau);

	// a = a1 + a2
	for(i = 0; i < dim; i++)
	{
		high_vec[i] = org_vec[i] + sigma; high_vec[i] -= sigma;
		low_vec[i] = org_vec[i] - high_vec[i];
	}
}

// error-free splitting
double ddotp_eft(double a[], double b[], int dim)
{
	int i;
	mpfr_t exact_val;
	double abs_a, abs_b, beta, v, w, tau, sigma;
	double a1[dim], a2[dim], b1[dim], b2[dim], a21[dim], a22[dim], b21[dim], b22[dim];
	double ip_a1b1, ip_a1b21, ip_a1b22, ip_a21b1, ip_a21b21, ip_a21b22, ip_a22b1, ip_a22b21, ip_a22b22;
	double ret, ip_array[4], ip_array_err[4];

	// (a1, a2) := split(a)
	// (b1, b2) := split(b)
	split_dvec(a1, a2, a, dim);
	split_dvec(b1, b2, b, dim);

	// (a21, a22) := split(a2)
	// (b21, b22) := split(b2)
	split_dvec(a21, a22, a2, dim);
	split_dvec(b21, b22, b2, dim);

	// ip_a1b1 = (a1, b1)
	// ip_a1b21 = (a1, b21)
	// ip_a, ip_a2b = (a2, b)
	ip_a1b1  = ddotp(a1, b1, dim);
	ip_a1b21 = ddotp(a1, b21, dim);
	ip_a1b22 = ddotp(a1, b22, dim);

	ip_a21b1  = ddotp(a21, b1, dim);
	ip_a21b21 = ddotp(a21, b21, dim);
	ip_a21b22 = ddotp(a21, b22, dim);

	ip_a22b1  = ddotp(a22, b1, dim);
	ip_a22b21 = ddotp(a22, b21, dim);
	ip_a22b22 = ddotp(a22, b22, dim);

/*	ip_array[0] = d_mpfrdotp_d(a1, b1, dim);
	ip_array[1] = d_mpfrdotp_d(a1, b2, dim);
	ip_array[2] = d_mpfrdotp_d(a2, b1, dim);
	ip_array[3] = d_mpfrdotp_d(a2, b2, dim);
*/

//	printf("a1b1, a1b2, a2b = %10.3e(%10.3e vs.%10.3e,(%10.3e, %10.3e), %10.3e, %10.3e\n", relerr_ddotp(a1, b1, dim), ip_a1b1, d_mpfrdotp_d(a1, b1, dim), ddotp(a1, a1, dim), ddotp(b1, b1, dim), relerr_ddotp(a1, b2, dim), relerr_ddotp(a2, b, dim));
//	printf("%10.3e(%10.3e vs.%10.3e), %10.3e, %10.3e, %10.3e\n", relerr_ddotp(a1, b1, dim), ip_a1b1, d_mpfrdotp_d(a1, b1, dim), relerr_ddotp(a2, b1, dim), relerr_ddotp(a2, b1, dim), relerr_ddotp(a2, b, dim));
//	printf(" %25.17e, %25.17e, %25.17e, %25.17e\n", ip_a1b1, ip_a1b2, ip_a2b1, ip_a2b2);

//	ret = qd::two_sum(ip_a1b2, ip_a2b, w);
//	ret = qd::two_sum(ip_a1b1, ret, v);
//	ret += v + w;

//	ret = ip_a1b1 + (ip_a1b21 + ip_a21b1) + (ip_a1b22 + ip_a21b21 + ip_a22b1) + (ip_a21b22 + ip_a22b21) + ip_a22b22;
//	ret = ip_a1b1 + (ip_a1b21 + ip_a21b1) + (ip_a1b22 + ip_a21b21 + ip_a22b1) + (ip_a21b22 + ip_a22b21) + ip_a22b22;
	ret = ip_a1b1;
	ip_array[0] = (ip_a1b21 + ip_a21b1);
	ip_array[1] = (ip_a1b22 + ip_a21b21 + ip_a22b1);
	ip_array[2] = (ip_a21b22 + ip_a22b21);
	ip_array[3] = ip_a22b22;

/*	ip_array[0] = qd::two_sum(ip_a1b21, ip_a21b1, ip_array_err[0]);
	ip_array[1] = qd::two_sum(ip_a1b22, ip_a21b21, ip_array_err[1]);
	ip_array[2] = qd::two_sum(ip_a1b22, ip_a21b21, ip_array_err[1]);
	ip_array[1] = (ip_a1b22 + ip_a21b21 + ip_a22b1);
	ip_array[2] = (ip_a21b22 + ip_a22b21);
	ip_array[3] = ip_a22b22;
*/
	qd::renorm(ret, ip_array[0], ip_array[1], ip_array[2], ip_array[3]);
/*	ret = ip_a1b1;
	ip_array[0] = qd::two_sum(ip_a1b21, ip_a21b1, ip_array_error[0]);
	ip_array[1] = qd::two_sum(ip_a1b22, ip_a21b21, ip_array_error[1]);
	ip_array[1] = qd::two_sum(ip_array[1], ip_a22b1, 
*/
//	ret = ip_a1b1 + ip_a1b2 + ip_a2b;
//	ret = ip_a1b1 + ip_a1b2 + ip_a2b1 + ip_a2b2;
//	ret = ip_a1b2 + ip_a2b1 + ip_a2b2;
//	ret = ip_a1b2 + ip_a2b;
//	ret = qd::two_sum(ip_a1b2, ip_a2b, v);
//	ret = qd::two_sum(ip_a2b1, ip_a2b2, v);
//	ret = qd::two_sum(ret, ip_a1b2 + v, w);
//	ret = d_mpfr_sum_d(ip_array, 4);
//	ret = qd::two_sum(ip_a1b2, ip_a2b, v);
//	ret = qd::two_sum(ip_a1b1, ret, w);
//	ret = ip_a1b1 + ip_a1b2;

	return ret;
}

// error-free splitting
double ddotp_eft2(double a[], double b[], int dim)
{
	int i;
	mpfr_t exact_val;
	double abs_a, abs_b, beta, v, w, tau, sigma;
	double a1[dim], a2[dim], b1[dim], b2[dim], a21[dim], a22[dim], a221[dim], a222[dim], b21[dim], b22[dim], b221[dim], b222[dim];
	double ip_a1b1, ip_a1b21, ip_a1b22, ip_a21b1, ip_a21b21, ip_a21b22, ip_a22b1, ip_a22b21, ip_a22b22;
	double ret, ip_array[4];

	// (a1, a2) := split(a)
	// (b1, b2) := split(b)
	split_dvec(a1, a2, a, dim);
	split_dvec(b1, b2, b, dim);

	// (a21, a22) := split(a2)
	// (b21, b22) := split(b2)
	split_dvec(a21, a22, a2, dim);
	split_dvec(b21, b22, b2, dim);

	// (a221, a222) := split(a22)
	// (b221, b222) := split(b22)
	split_dvec(a221, a222, a22, dim);
	split_dvec(b221, b222, b22, dim);

	// ip_a1b1 = (a1, b1)
	// ip_a1b21 = (a1, b21)
	// ip_a, ip_a2b = (a2, b)
	ip_a1b1  = ddotp(a1, b1, dim);
	ip_a1b21 = ddotp(a1, b21, dim);
	//ip_a1b22 = ddotp(a1, b22, dim);
	ip_a1b22 = ddotp(a1, b221, dim) + ddotp(a1, b222, dim);

	ip_a21b1  = ddotp(a21, b1, dim);
	ip_a21b21 = ddotp(a21, b21, dim);
//	ip_a21b22 = ddotp(a21, b22, dim);
	ip_a21b22 = ddotp(a21, b221, dim) + ddotp(a21, b222, dim);

//	ip_a22b1  = ddotp(a22, b1, dim);
	ip_a22b1  = ddotp(a221, b1, dim) + ddotp(a222, b1, dim);
//	ip_a22b21 = ddotp(a22, b21, dim);
	ip_a22b21 = ddotp(a221, b21, dim) + ddotp(a222, b21, dim);
//	ip_a22b22 = ddotp(a22, b22, dim);
	ip_a22b22 = ddotp(a221, b221, dim) + ddotp(a221, b222, dim) + ddotp(a221, b221, dim) + ddotp(a222, b222, dim);

/*	ip_array[0] = d_mpfrdotp_d(a1, b1, dim);
	ip_array[1] = d_mpfrdotp_d(a1, b2, dim);
	ip_array[2] = d_mpfrdotp_d(a2, b1, dim);
	ip_array[3] = d_mpfrdotp_d(a2, b2, dim);
*/

//	printf("a1b1, a1b2, a2b = %10.3e(%10.3e vs.%10.3e,(%10.3e, %10.3e), %10.3e, %10.3e\n", relerr_ddotp(a1, b1, dim), ip_a1b1, d_mpfrdotp_d(a1, b1, dim), ddotp(a1, a1, dim), ddotp(b1, b1, dim), relerr_ddotp(a1, b2, dim), relerr_ddotp(a2, b, dim));
//	printf("%10.3e(%10.3e vs.%10.3e), %10.3e, %10.3e, %10.3e\n", relerr_ddotp(a1, b1, dim), ip_a1b1, d_mpfrdotp_d(a1, b1, dim), relerr_ddotp(a2, b1, dim), relerr_ddotp(a2, b1, dim), relerr_ddotp(a2, b, dim));
//	printf(" %25.17e, %25.17e, %25.17e, %25.17e\n", ip_a1b1, ip_a1b2, ip_a2b1, ip_a2b2);

//	ret = qd::two_sum(ip_a1b2, ip_a2b, w);
//	ret = qd::two_sum(ip_a1b1, ret, v);
//	ret += v + w;

//	ret = ip_a1b1 + (ip_a1b21 + ip_a21b1) + (ip_a1b22 + ip_a21b21 + ip_a22b1) + (ip_a21b22 + ip_a22b21) + ip_a22b22;
	ret = ip_a1b1 + (ip_a1b21 + ip_a21b1) + (ip_a1b22 + ip_a21b21 + ip_a22b1) + (ip_a21b22 + ip_a22b21) + ip_a22b22;
/*	ret = ip_a1b1;
	ip_array[0] = qd::two_sum(ip_a1b21, ip_a21b1, ip_array_error[0]);
	ip_array[1] = qd::two_sum(ip_a1b22, ip_a21b21, ip_array_error[1]);
	ip_array[1] = qd::two_sum(ip_array[1], ip_a22b1, 
*/
//	ret = ip_a1b1 + ip_a1b2 + ip_a2b;
//	ret = ip_a1b1 + ip_a1b2 + ip_a2b1 + ip_a2b2;
//	ret = ip_a1b2 + ip_a2b1 + ip_a2b2;
//	ret = ip_a1b2 + ip_a2b;
//	ret = qd::two_sum(ip_a1b2, ip_a2b, v);
//	ret = qd::two_sum(ip_a2b1, ip_a2b2, v);
//	ret = qd::two_sum(ret, ip_a1b2 + v, w);
//	ret = d_mpfr_sum_d(ip_array, 4);
//	ret = qd::two_sum(ip_a1b2, ip_a2b, v);
//	ret = qd::two_sum(ip_a1b1, ret, w);
//	ret = ip_a1b1 + ip_a1b2;

	return ret;
}


double drand(void)
{
	return (double)rand() / RAND_MAX;
}

// GenDot
// ref. Ogita et.al., "Accurate sum and dot product", SIAM 
double dgendot(double x[], double y[], double &true_cond, int n, double cond)
{
	int i, n2 = n / 2, index;
	double b, exponent[n], h, true_dotp;
	mpfr_t mpfr_exact_dotp;
	mpfr_prec_t exact_dotp_prec;

	// initialize
	b = log2(cond);
//	printf("b = %25.17e\n", b);

//	mpfr_init2(mpfr_exact_dotp, (mpfr_prec_t)ceil(cond * 2));

	exact_dotp_prec = (mpfr_prec_t)ceil(b * 1.5);
	if(exact_dotp_prec < 64) exact_dotp_prec = 64;

//	printf("cond, prec = %10.3e, %ld\n", cond, exact_dotp_prec);

	mpfr_init2(mpfr_exact_dotp, exact_dotp_prec);

	for(i = 0; i < n; i++)
	{
		x[i] = 0.0;
		y[i] = 0.0;
	}

	for(i = 0; i < n2; i++)
		exponent[i] = round(drand() * b / 2);

	exponent[0] = round(b / 2) + 1;
	exponent[n2 - 1] = 0;

	h = (b / 2) / (double)(n - n2);
	for(i = 0; i < n2; i++)
	{
		x[i] = (2 * drand() - 1) * pow(2.0, exponent[i]);
		y[i] = (2 * drand() - 1) * pow(2.0, exponent[i]);

		exponent[i + n2] = round(h * (i + 1));
	}
//	printf("start x, y set up!\n");

	for(i = n2 ; i < n; i++)
	{
		x[i] = (2 * drand() - 1) * pow(2.0, exponent[i]);
		mpfr_mydotp_d(mpfr_exact_dotp, x, y, i);
		true_dotp = mpfr_get_d(mpfr_exact_dotp, MPFR_RNDN);
		y[i] = ((2 * drand() - 1) * pow(2.0, exponent[i]) - true_dotp) / x[i];
	}

//	for(i = 0; i < n; i++)
//		printf("exponent[%d] = %10.3e, x[%d] = %25.17e\n", i, exponent[i], i, x[i]);

//	printf("x, y have been set up!\n");

	// random permulation
	for(i = 0; i < n2; i++)
	{
		index = rand() % n;

		//swap(x[rand() % n], x[i])
		h = x[i];
		x[i] = x[index];
		x[index] = h;

		h = y[i];
		y[i] = y[index];
		y[index] = h;
	}

	mpfr_mydotp_d(mpfr_exact_dotp, x, y, n);
	true_dotp = mpfr_get_d(mpfr_exact_dotp, MPFR_RNDN);

//	printf("true_dot_p = %25.17e\n", true_dotp);

	mpfr_clear(mpfr_exact_dotp);

	// true condition number
	true_cond = 0.0;
	for(i = 0; i < n; i++)
		true_cond += fabs(x[i]) * fabs(y[i]);

	true_cond = 2 * true_cond / fabs(true_dotp);

	return true_dotp;
}

// simple gendot
// (x, y) := ||x||_2 * ||y||_2 * cos(theta)
// -> x, y = random()/||x||_2 or ||y||_2 -> y' := y - (x, y) * x / ||x|| * cos(theta)
double simple_dgendot(double x[], double y[], double &true_cond, int n, double cond)
{
	int i, n2 = n/2;
	double b, exponent[n], nx[n], h, true_dotp, norm2_x, ip_xy, sin_theta, theta;
	mpfr_t mpfr_exact_dotp;
	mpfr_prec_t exact_dotp_prec;

	// initialize
	b = log2(cond);
//	theta = acos(1/fabs(cond)); // true cond > cond
	theta = acos(1/sqrt(fabs(cond))); // true_cond < cond
//	printf("b = %25.17e\n", b);

//	mpfr_init2(mpfr_exact_dotp, (mpfr_prec_t)ceil(cond * 2));

	exact_dotp_prec = (mpfr_prec_t)ceil(b * 1.5);
	if(exact_dotp_prec < 64) exact_dotp_prec = 64;

//	printf("cond, cos(theta), prec = %10.3e, %10.3e, %ld\n", cond, cos(theta), exact_dotp_prec);

	mpfr_init2(mpfr_exact_dotp, exact_dotp_prec);

	for(i = 0; i < n; i++)
		exponent[i] = round(drand() * b / 2);

	for(i = 0; i < n; i++)
	{
		x[i] = (2 * drand() - 1) * pow(2.0, exponent[i]);
		y[i] = (2 * drand() - 1) * pow(2.0, exponent[i]);
	}

	// normalize
	norm2_x = ddotp(x, x, n);
	norm2_x = sqrt(norm2_x);

	for(i = 0 ; i < n; i++)
		nx[i] = x[i] / norm2_x;

	ip_xy = ddotp(nx, y, n);
	sin_theta = sin(theta);

	for(i = 0; i < n; i++)
		y[i] -= ip_xy * nx[i] * sin_theta;

//	for(i = 0; i < n; i++)
//		printf("exponent[%d] = %10.3e, x[%d] = %25.17e\n", i, exponent[i], i, x[i]);

//	printf("x, y have been set up!\n");

	mpfr_mydotp_d(mpfr_exact_dotp, x, y, n);
	true_dotp = mpfr_get_d(mpfr_exact_dotp, MPFR_RNDN);

//	printf("true_dot_p = %25.17e\n", true_dotp);

	mpfr_clear(mpfr_exact_dotp);

	// true condition number
	true_cond = 0.0;
	for(i = 0; i < n; i++)
		true_cond += fabs(x[i]) * fabs(y[i]);

	true_cond = 2 * true_cond / fabs(true_dotp);

	return true_dotp;
}


int main()
{
	int dim;
	double init_cond;
	double a, b, c, err;
	double a_hi, a_low;

	double *x, *y, cond, dotp;

// double

//	a = sqrt(2.0);
//	a = -sqrt(2.0) + sqrt(3.0);
//	a = pow(2.0, (double)53) - 1; // full 53-bit ones
	a = pow(2.0, 53.0) -  pow(2.0, 26.0) - 1; // full 53-bit ones
//	a = pow(2.0, 53.0) -  pow(2.0, 25.0) - 1; // full 53-bit ones
//	a = pow(2.0, 53.0) -  pow(2.0, 27.0) - 1; // full 53-bit ones
//	a = pow(2.0, 53.0) -  pow(2.0, 28.0) - 1; // full 53-bit ones
//	a = pow(2.0, 53.0) -  pow(2.0, 29.0) - 1; // full 53-bit ones
	b = sqrt(3.0);

	// quick_two_sum
	c = dquick_two_sum(a, b, err);
	printf("QuickTwoSum(%25.17e, %25.17e) = %25.17e, %25.17e\n", a, b, c, err);
	printf("a = "); printb_double(a);
	printf("b = "); printb_double(b);
	printf("c = "); printb_double(c);
	printf("err = "); printb_double(err);

	// split
	split_d2d(a_hi, a_low, a);
	printf("Split(%25.17e) = %25.17e, %25.17e\n", a, a_hi, a_low);
	printb_double(a);
	printb_double(a_hi); printb_double(a_low);

// float

	float af, bf, cf, errf;
	float af_high, af_low;

	af = powf(2.0, 24.0) -  powf(2.0, 13.0) - 1; // full 53-bit ones
	bf = sqrtf(3.0);

	// quick_two_sum
	cf = fquick_two_sum(af, bf, errf);
	printf("QuickTwoSum(%15.7e, %15.7e) = %15.7e, %15.7e\n", af, bf, cf, errf);
	printf("af = "); printb_float(af);
	printf("bf = "); printb_float(bf);
	printf("cf = "); printb_float(cf);
	printf("errf = "); printb_float(errf);

	// split
	split_f2f(af_high, af_low, af);
	printf("Split(%15.7e) = %15.7e, %15.7e\n", af, af_high, af_low);
	printb_float(af);
	printb_float(af_high); printb_float(af_low);

	// generate ill-condtioned vectors
	for(dim = 10; dim <= 10000; dim *= 10)
//	for(dim = 1; dim <= 2; dim++)
	{
		x = (double *)calloc(dim, sizeof(double));
		y = (double *)calloc(dim, sizeof(double));

		printf("----- dim = %d ----- \n", dim);
		printf("dotp      , init_cond : cond  , ddotp_relerr  dotpxblas_error\n");
		for(init_cond = 10.0; init_cond < 1.0e+30; init_cond *= 100)
		{
			dotp = dgendot(x, y, cond, dim, init_cond);
			//dotp = simple_dgendot(x, y, cond, dim, init_cond);
			printf("%10.3e : %10.3e, %10.3e, %10.3e, %10.3e\n", init_cond, cond, \
				fabs((ddotp(x, y, dim) - dotp) / dotp), \
				fabs((ddotp_xblas(x, y, dim) - dotp) / dotp), \
				fabs((ddotp_eft(x, y, dim) - dotp) / dotp)
//				fabs((ddotp_eft2(x, y, dim) - dotp) / dotp)
			);
		
		//	printf("%25.17e, %25.17e, %25.17e\n", ddotp(x, y, dim), ddotp_xblas(x, y, dim), ddotp_eft(x, y, dim));
		}
		for(int i = 0; i< dim; i++)
		{
			//printf("%5d, %25.17e, %25.17e\n", i, x[i], y[i]);
		//	printb_double(x[i]); printb_double(y[i]);
		}

		free(x);
		free(y);
	}

//	for(int i = 0; i < 10; i++)
//		printf("x[%d], y[%d] = %25.17e, %25.17e\n", i, i, x[i], y[i]);

	return 0;
}
