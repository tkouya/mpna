#include <iostream>
#include <cmath>
#include <cstdio>
#include <qd/inline.h>

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

int main()
{
	double a, b, c, err;
	double a_hi, a_low;

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

//	af = powf(2.0, 24.0) -  powf(2.0, 13.0) - 1; // full 53-bit ones
//	af = powf(2.0, 24.0) -  powf(2.0, 12.0) - 1; // full 53-bit ones
	af = powf(2.0, 24.0) -  powf(2.0, 11.0) - 1; // full 53-bit ones
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

	return 0;
}
