// logistic ŽÊ‘œ
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// true values by logistic_qd.cc
double true_x[] = {
	7.50100000000000000e-01, // x[ 0]
	7.49799960000000000e-01, // x[ 1]
	7.50399919935993600e-01, // x[ 2]
	7.49199520384191979e-01, // x[ 3]
	7.51598396161154744e-01, // x[ 4]
	7.46792988196538534e-01, // x[ 5]
	7.56372883908092767e-01, // x[ 6]
	7.37091777586590356e-01, // x[ 7]
	7.75149956003323084e-01, // x[ 8]
	6.97170006845477485e-01, // x[ 9]
	8.44495953602217448e-01, // x[10]
	5.25290151806795375e-01, // x[11]
	9.97441632886356979e-01, // x[12]
	1.02072874854194051e-02, // x[13]
	4.04123950704376218e-02, // x[14]
	1.55116933580433963e-01, // x[15]
	5.24222681988148806e-01, // x[16]
	9.97653046709204046e-01, // x[17]
	9.36578040418710571e-03, // x[18]
	3.71122502464306021e-02, // x[19]
	1.42939724512307655e-01, // x[20]
	4.90031838674612997e-01, // x[21]
	9.99602543039164235e-01, // x[22]
	1.58919595520019145e-03, // x[23]
	6.34668164566466722e-03, // x[24]
	2.52256051110130018e-02, // x[25]
	9.83570958311849461e-02, // x[26]
	3.54731910123360188e-01, // x[27]
	9.15588728254369991e-01, // x[28]
	3.09144035791661656e-01, // x[29]
	8.54296003704421892e-01, // x[30]
	4.97897367036305077e-01, // x[31]
	9.99982315738479934e-01, // x[32]
	7.07357951478418944e-05, // x[33]
	2.82923166380506788e-04, // x[34]
	1.13137248344972806e-03, // x[35]
	4.52036991901368344e-03, // x[36]
	1.79997446992358386e-02, // x[37]
	7.07030155599926799e-02, // x[38]
	2.62816396602864452e-01, // x[39]
	7.74975753118201241e-01, // x[40]
	6.97553340788312161e-01, // x[41]
	8.43890710173507994e-01, // x[42]
	5.26956717825441302e-01, // x[43]
	9.97093341456318141e-01, // x[44]
	1.15928395191692022e-02, // x[45]
	4.58337823642079640e-02, // x[46]
	1.74932187033593532e-01, // x[47]
	5.77323667892949532e-01, // x[48]
	9.76084201534323383e-01, // x[49]
	9.33753321977030291e-02, // x[50]
	3.38625518138686532e-01, // x[51]
	8.95833106415970442e-01, // x[52]
	3.73264607460332078e-01, // x[53]
	9.35752561111265137e-01, // x[54]
	2.40478821939892567e-01, // x[55]
	7.30595032553176050e-01, // x[56]
	7.87303723847198710e-01, // x[57]
	6.69826281054130332e-01, // x[58]
	8.84636137053294132e-01, // x[59]
	4.08220168290878132e-01, // x[60]
	9.66305849965781072e-01, // x[61]
	1.30235417150761889e-01, // x[62]
	4.53096613081315701e-01, // x[63]
	9.91200289182224778e-01, // x[64]
	3.48891036291950051e-02, // x[65]
	1.34687416308585188e-01, // x[66]
	4.66186864786812195e-01, // x[67]
	9.95426687548218716e-01, // x[68]
	1.82095890599986669e-02, // x[69]
	7.15119997050585749e-02, // x[70]
	2.65592134412969109e-01, // x[71]
	7.80211810203729835e-01, // x[72]
	6.85925365689395552e-01, // x[73]
	8.61727033573058140e-01, // x[74]
	4.76614212729742676e-01, // x[75]
	9.97812419815001082e-01, // x[76]
	8.73117871173247273e-03, // x[77]
	3.46197809201450496e-02, // x[78]
	1.33685006756744842e-01, // x[79]
	4.63253302900775711e-01, // x[80]
	9.94598721009191445e-01, // x[81]
	2.14884207042880214e-02, // x[82]
	8.41066739196941891e-02, // x[83]
	3.08130965287441690e-01, // x[84]
	8.52745094073884377e-01, // x[85]
	5.02283594427225845e-01, // x[86]
	9.99979140785967772e-01, // x[87]
	8.34351157016707688e-05, // x[88]
	3.33712617132554471e-04, // x[89]
	1.33440501208688405e-03, // x[90]
	5.33049750140240580e-03, // x[91]
	2.12083331911597940e-02, // x[92]
	8.30341591776501754e-02, // x[93]
	3.04557950349243314e-01, // x[94]
	8.47209620913244634e-01, // x[95]
	5.17781916581123819e-01, // x[96]
	9.98735213770807815e-01, // x[97]
	5.05274617994652239e-03, // x[98]
	2.01088637439502329e-02, // x[99]
	7.88179893715099068e-02  // x[100]
};

double relerr(double approx, double true_val)
{
	double relative_error = approx - true_val;

	if(true_val != 0)
		relative_error /= true_val;
	
	relative_error = fabs(relative_error);

	return relative_error;
}

// double prec. 
double dlogistic(double x_init, int print_interval, int max_k)
{
	int k;
	double x_kp1, x_k;

	// x0 := x_init
	x_k = x_init;

	for(k = 0; k < max_k; k++)
	{
		// x_{k+1} := 4 * x_k * (1 - x_k)
		x_kp1 = 4 * x_k * (1 - x_k);

		if((k % print_interval) == 0)
			printf("%5d, %25.17e, %10.3e\n", k, x_k, relerr(x_k, true_x[k]));

		// x_k := x_{k+1}
		x_k = x_kp1;
	}

	return x_k;
}

// 2-fold precision

#ifndef DDSIZE
	#define DDSIZE 2
#endif // DDSIZE

#ifndef QDSIZE
	#define QDSIZE 4
#endif // QDSIZE

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
double quick_two_sum(double a, double b, double *err)
{
	double s = a + b;
	*err = b - (s - a);
	return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
double quick_two_diff(double a, double b, double *err)
{
	double s = a - b;
	*err = (a - s) - b;
	return s;
}

/* Computes fl(a+b) and err(a+b).  */
double two_sum(double a, double b, double *err)
{
	double s = a + b;
	double bb = s - a;
	*err = (a - (s - bb)) + (b - bb);
	return s;
}

/* Computes fl(a-b) and err(a-b).  */
double two_diff(double a, double b, double *err)
{
	double s = a - b;
	double bb = s - a;
	*err = (a - (s - bb)) - (b + bb);
	return s;
}

#ifndef QD_FMS
/* Computes high word and lo word of a */
void split(double a, double *hi, double *lo)
{
	double temp;

	if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH)
	{
		a *= 3.7252902984619140625e-09;  // 2^-28
		temp = _QD_SPLITTER * a;
		*hi = temp - (temp - a);
		*lo = a - (*hi);
		*hi *= 268435456.0;          // 2^28
		*lo *= 268435456.0;          // 2^28
	}
	else
	{
		temp = _QD_SPLITTER * a;
		*hi = temp - (temp - a);
		*lo = a - (*hi);
	}
}
#endif // QD_FMS

// check FP_FAST_FMA
#ifdef FP_FAST_FMA
	#define QD_FMA(a, b, c) fma((a), (b), (c)) // a * b + c
	#define QD_FMS(a, b, c) QD_FMA((a), (b), -(c)) // a * b - c
#endif // FP_FAST_FMA


/* Computes fl(a*b) and err(a*b). */
double two_prod(double a, double b, double *err)
{
#ifdef QD_FMS
	double p = a * b;

	*err = QD_FMS(a, b, p);

	return p;
#else
	double a_hi, a_lo, b_hi, b_lo;
	double p = a * b;

	split(a, &a_hi, &a_lo);
	split(b, &b_hi, &b_lo);
	*err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;

	return p;
#endif
}

#ifdef X86
#ifdef  _WIN32
#include <float.h>
#else

#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif

#endif
#endif /* X86 */

void fpu_fix_start(unsigned int *old_cw) {
#ifdef X86
#ifdef _WIN32
#ifdef __BORLANDC__
  /* Win 32 Borland C */
  unsigned short cw = _control87(0, 0);
  _control87(0x0200, 0x0300);
  if (old_cw) {
    *old_cw = cw;
  }
#else
  /* Win 32 MSVC */
  unsigned int cw = _control87(0, 0);
  _control87(0x00010000, 0x00030000);
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#else
  /* Linux */
  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);

  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
  
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#endif
}

void fpu_fix_end(unsigned int *old_cw) {
#ifdef X86
#ifdef _WIN32

#ifdef __BORLANDC__
  /* Win 32 Borland C */
  if (old_cw) {
    unsigned short cw = (unsigned short) *old_cw;
    _control87(cw, 0xFFFF);
  }
#else
  /* Win 32 MSVC */
  if (old_cw) {
    _control87(*old_cw, 0xFFFFFFFF);
  }
#endif

#else
  /* Linux */
  if (old_cw) {
    int cw;
    cw = *old_cw;
    _FPU_SETCW(cw);
  }
#endif
#endif
}
/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
double two_sqr(double a, double *err)
{
#ifdef QD_FMS
	double p = a * a;

	*err = QD_FMS(a, a, p);

	return p;
#else
	double hi, lo;
	double q = a * a;

	split(a, &hi, &lo);
	*err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;

	return q;
#endif
}

// double prec. 
/* [tkouya@hpcsv01 mpna]$ icc -fp-model precise -DFP_FAST_FMA logistic_2f.c -lm
[tkouya@hpcsv01 mpna]$ ./a.out
    0,   7.50099999999999989e-01,  0.000e+00
   10,   8.44495953602201199e-01,  1.919e-14
   20,   1.42939724528399537e-01,  1.126e-10
   30,   8.54296020314658677e-01,  1.944e-08
   40,   7.74995885155205677e-01,  2.598e-05
   50,   7.95128764501052410e-02,  1.485e-01
   60,   2.73187240440892098e-01,  3.308e-01
   70,   5.52530562083362264e-01,  6.726e+00
   80,   2.16255663995813446e-01,  5.332e-01
   90,   7.87467937188412348e-01,  5.891e+02
  100,   2.69706745887651977e-01,  2.422e+00

    0,   7.50099999999999989e-01 ( 0.000e+00)
   10,   8.44495953602217408e-01 ( 0.000e+00)
   20,   1.42939724512307659e-01 ( 0.000e+00)
   30,   8.54296003704421936e-01 ( 0.000e+00)
   40,   7.74975753118201216e-01 ( 0.000e+00)
   50,   9.33753321977030259e-02 ( 0.000e+00)
   60,   4.08220168290875429e-01 ( 6.663e-15)
   70,   7.15119997065230134e-02 ( 2.048e-11)
   80,   4.63253305802721616e-01 ( 6.264e-09)
   90,   1.33462256598348853e-03 ( 1.630e-04)
  100,   7.71815430797299223e-02 ( 2.076e-02)
*/
double dlogistic_2fold(double x_init, double x_init_err, int print_interval, int max_k)
{
	int k;
	double x_kp1, x_k, w[3], t, e_kp1, e_k, e_w[3];

	// x0 := x_init
	x_k = x_init;
	
	// e_0 := 0;
	//e_k = 0.0;
	e_k = x_init_err;

	for(k = 0; k < max_k; k++)
	{
		// x_k := x_k + e_k
		//x_k += e_k;

		// x_{k+1} := 4 * x_k * (1 - x_k)

		// (w, e_w) :=  (1 - x_k);
		w[0] = two_diff(1.0, x_k, &e_w[0]);
		//w = 1 - x_k;

		// e_k := e_w - e_k
		//e_k = e_w[0] - e_k;
		//printf("e_k, e_w = %25.17e, %25.17e\n", e_k, e_w);

		// (w, e_w) := x_k * (1 - x_k)
		w[1] = two_prod(x_k, w[0], &e_w[1]);
		//w = x_k * w;

		// e_k := e_w + x_k * e_k
		//e_k = e_w[1] + x_k * e_k;
		//printf("e_k, e_w = %25.17e, %25.17e\n", e_k, e_w);

		// (w, e_w) := 4 * x_k * (1 - x_k)
		w[2] = two_prod(4.0, w[1], &e_w[2]);
		//w = 4.0 * w;
	
		// e_k := e_w + 4 * e_k
		//e_k = e_w[2] + 4 * e_k;
		//printf("e_k, e_w = %25.17e, %25.17e\n", e_k, e_w);

		// e_{k+1} := 4 * x_k * e_1 + 4 * e_2 + e_3 - 4 * x_k * e_k + 4 * (1 - x_k) * e_k - 4 * e_k^2
		//e_kp1 = 4 * x_k * e_w[0] + 4 * e_w[1] + e_w[2] - 4 * x_k * e_k + 4 * w[0] * e_k - 4 * e_k * e_k;
		e_kp1 = 4 * x_k * e_w[0] + 4 * e_w[1] + e_w[2] - 4 * x_k * e_k + 4 * w[0] * e_k;
		

		// (x_{k+1}, e_w) := x_{k+1} + e_k
		//x_kp1 = two_sum(w[2], e_kp1, &e_w[0]);
		x_kp1 = quick_two_sum(w[2], e_kp1, &e_w[0]);
		//x_kp1 = quick_two_sum(w[2], e_k, &e_w[0]);
		//x_kp1 = w + e_k;
		//x_kp1 = w;
		//printf("e_k, e_w = %25.17e, %25.17e\n", e_k, e_w);

		if( (k % print_interval) == 0)
			printf("%5d, %25.17e (%10.3e)\n", k, x_k, relerr(x_k, true_x[k]));

		// x_k := x_{k+1}
		x_k = x_kp1;
		e_k = e_w[0];

	}

	return x_k;
}



int main()
{
	int i;
	double x0, x0_err, x_final;

	// ‰Šú’l
	//x0 = 0.7501;
	//x0     = 0.75009999999999999;
	//x0_err = 1.1013412404281551e-17;
	x0     = 0.75009999999999998899;
	x0_err = 1.1013412404281550953e-17;
	//x0_err = 0.0;

	// x_100
	x_final = dlogistic(x0, 10, 101);

	fpu_fix_start(NULL);

	x_final = dlogistic_2fold(x0, x0_err, 10, 101);

	return 0;
}
