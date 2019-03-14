#include <stdio.h>
#include <stdlib.h>

#include "mpfr.h"

// func(x) = 0

// func(x) 
void mpfr_func(mpfr_t ret, mpfr_t x, mpfr_rnd_t rmode)
{
	// ret := x - cos(x)
	mpfr_cos(ret, x, rmode); mpfr_sub(ret, x, ret, rmode);

	// ret := x^2 - a
	//mpfr_sqr(ret, x, rmode); mpfr_sub_ui(ret, ret, 2, rmode);
}

// dfunc(x) = func(x)' = 
void mpfr_dfunc(mpfr_t ret, mpfr_t x, mpfr_rnd_t rmode)
{
	// ret := 1 + sin(x)
	mpfr_sin(ret, x, rmode); mpfr_add_ui(ret, ret, 1UL, rmode);

	// ret := 2 * x
	//mpfr_mul_ui(ret, x, 2UL, rmode);
}

int main(int argc, char *argv[])
{
	mp_bitcnt_t old_prec, new_prec, max_prec;
	mpfr_t old_x, new_x, func_x, dfunc_x, relerr;
	int itimes;

	if(argc <= 1)
	{
		printf("Usage: %s [max_prec in bits] \n", argv[0]);
		return 0;
	}

	max_prec = (mp_bitcnt_t)atoi(argv[1]);
	if(max_prec <= 2)
	{
		fprintf(stderr, "Error: max_prec = %ld is too small!\n", max_prec);
		return -1;
	}

	old_prec = 2; // 1bit

	mpfr_init2(old_x, old_prec);
	mpfr_init2(new_x, old_prec);
	mpfr_init2(func_x, old_prec);
	mpfr_init2(dfunc_x, old_prec);
	mpfr_init2(relerr, max_prec);

	mpfr_set_ui(old_x, 1UL, MPFR_RNDN);
	//mpfr_set_d(old_x, 0.9, MPFR_RNDN);
	//mpfr_set_d(old_x, 0.3, MPFR_RNDN);

	// newton iteration
	for(itimes = 1; itimes <= 20; itimes++)
	{
		// new_x := old_x - func(old_x) / dfunc(old_x)
		mpfr_func(func_x, old_x, MPFR_RNDN); mpfr_dfunc(dfunc_x, old_x, MPFR_RNDN); mpfr_div(new_x, func_x, dfunc_x, MPFR_RNDN);
		mpfr_sub(new_x, old_x, new_x, MPFR_RNDN);
		
		// new_x := x * (2 - a * x)
		//mpfr_mul_ui(new_x, old_x, 3UL, MPFR_RNDN); mpfr_ui_sub(new_x, 2UL, new_x, MPFR_RNDN); mpfr_mul(new_x, old_x, new_x, MPFR_RNDN);

		// detect precision
		
		mpfr_reldiff(relerr, new_x, old_x, MPFR_RNDN);
		//mpfr_printf("relerr = %RNe\n", relerr);
/*		if(mpfr_cmp_ui(relerr, 0UL) > 0)
		{
			mpfr_log10(relerr, relerr, MPFR_RNDN);
			new_prec = (mp_bitcnt_t)ceil(-mpfr_get_d(relerr, MPFR_RNDN) / log10(2.0));
			mpfr_printf("new_prec = %ld, relerr = %RNe\n", new_prec, relerr);
		}
		else
			new_prec = max_prec;
*/
		new_prec = 2 * old_prec;

		mpfr_printf("%3d %10.3RNe %50.43RNe\n", itimes, relerr, new_x);

		// reset prec
//		if(old_prec >= max_prec)
//			break;

		printf("old_prec -> new_prec : max_prec= %ld -> %ld: %ld\n", old_prec, new_prec, max_prec);

		mpfr_prec_round(new_x, new_prec, MPFR_RNDN);
		mpfr_prec_round(old_x, new_prec, MPFR_RNDN);
		mpfr_prec_round(func_x, new_prec, MPFR_RNDN);
		mpfr_prec_round(dfunc_x, new_prec, MPFR_RNDN);
		old_prec = new_prec;

		// old_x := new_x
		mpfr_set(old_x, new_x, MPFR_RNDN);
	}

	mpfr_clear(old_x);
	mpfr_clear(new_x);
	mpfr_clear(func_x);
	mpfr_clear(dfunc_x);
	mpfr_clear(relerr);

	return 0;
}
