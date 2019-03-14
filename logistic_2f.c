// logistic ŽÊ‘œ
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "mpfr.h"

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
			printf("%5d, %25.17e, %10.3e\n", k, x_k);
			//printf("%5d, %25.17e, %10.3e\n", k, x_k, relerr(x_k, true_x[k]));

		// x_k := x_{k+1}
		x_k = x_kp1;
	}

	return x_k;
}


// mpfr 
void mpfr_logistic(mpfr_t ret, const char *x_init_str, int print_interval, int max_k, unsigned long prec)
{
	int k;
	mpfr_t x_kp1, x_k;

	mpfr_init2(x_kp1, prec);
	mpfr_init2(x_k, prec);

	// x0 := x_init
	//x_k = x_init;
	mpfr_set_str(x_k, x_init_str, 10, MPFR_RNDN);

	for(k = 0; k < max_k; k++)
	{
		// x_{k+1} := 4 * x_k * (1 - x_k)
		//x_kp1 = 4 * x_k * (1 - x_k);
		mpfr_ui_sub(x_k, 1UL, x_k, MPFR_RNDN);
		mpfr_mul(x_kp1, x_kp1, x_k, MPFR_RNDN);
		mpfr_mul_ui(x_kp1, x_kp1, 4UL, MPFR_RNDN);

		//if((k % print_interval) == 0)
			//mpfr_printf("%5d, %25.17RNe\n", k, x_k);
			mpfr_printf("%5d, %25.17RNe %25.17RNe\n", k, x_k, x_kp1);

		// x_k := x_{k+1}
		mpfr_set(x_k, x_kp1, MPFR_RNDN);
	}

	mpfr_set(ret, x_kp1, MPFR_RNDN);

	mpfr_clear(x_kp1);
	mpfr_clear(x_k);

	return;
}

double relerr_mpfr(double approx, mpfr_t true_val)
{
	mpfr_t tmp, ret;

	mpfr_init2(tmp, mpfr_get_prec(true_val));
	mpfr_init2(ret, mpfr_get_prec(true_val));

	mpfr_set_d(tmp, approx, MPFR_RNDN);

	mpfr_reldiff(ret, true_val, tmp, MPFR_RNDN);

	mpfr_clear(tmp);
	mpfr_clear(ret);
}

// 2-fold 
#include "etf.c"

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

	mpfr_t true_val;
	unsigned long prec = 256;

	mpfr_init2(true_val, prec);

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
		{
			//printf("%5d, %25.17e (%10.3e)\n", k, x_k, relerr(x_k, true_x[k]));
			mpfr_logistic(true_val, "0.7501", 0, k, prec);
			//printf("%5d, %25.17e (%10.3e)\n", k, x_k, relerr_mpfr(x_k, true_val));
			mpfr_printf("%5d, %25.17e %25.17RNe (%10.3e)\n", k, x_k, true_val, relerr_mpfr(x_k, true_val));
		}

		// x_k := x_{k+1}
		x_k = x_kp1;
		e_k = e_w[0];

	}

	mpfr_clear(true_val);

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
