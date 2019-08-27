#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpfr.h"

// Relative error
void relerr(mpfr_t c, mpfr_t a, mpfr_t ad)
{
	mpfr_sub(c, a, ad, MPFR_RNDN);
	if(mpfr_cmp_ui(ad, 0UL) != 0)
		mpfr_div(c, c, ad, MPFR_RNDN);
	mpfr_abs(c, c, MPFR_RNDN);
}

// mpfr_myexp(x)
void mpfr_myexp(mpfr_t ret, mpfr_t x, int n, mpfr_rnd_t rmode)
{
	int i;
	mpfr_prec_t prec;
	mpfr_t *coef, xn;

	// set default prec in this function
	prec = mpfr_get_prec(ret);

	// initialize
	mpfr_init2(xn, prec);

	coef = (mpfr_t *)calloc(n + 1, sizeof(mpfr_t));
	for(i = 0; i <= n; i++)
		mpfr_init2(coef[i], prec);

	// coef[i] = 1 / i!
	mpfr_set_ui(coef[0], 1UL, rmode);
	for(i = 1; i <= n; i++)
		mpfr_div_ui(coef[i], coef[i - 1], i, rmode);

	// Horner method
	mpfr_set(ret, coef[n], rmode);
	for(i = n - 1; i >= 0; i--)
		mpfr_fma(ret, ret, x, coef[i], rmode);

	// clear
	mpfr_clear(xn);
	for(i = 0; i <= n; i++)
		mpfr_clear(coef[i]);
	free(coef);
}

int main()
{
	int n, str_x_length, i;
	char *str_x;
	mpfr_prec_t prec, long_prec;
	mpfr_t mp_myexp, mp_myexp_long, mp_x, mp_exp_long, mp_myexp_relerr, mp_myexp_long_relerr, mp_myexp_rf_relerr;

	printf("Input prec in bits: "); while(scanf("%ld", &prec) < 1);
	// printf("prec  = %ld\n", prec);
	long_prec = (prec < 64) ? 256 : prec * 3;

	printf("#terms : "); while(scanf("%d", &n) < 1);
	// printf("#term = %d\n", n);

	str_x_length = (int)ceil(log10(2.0) * prec);
	str_x = (char *)calloc(str_x_length, sizeof(char));

	mpfr_init2(mp_myexp, prec);
	mpfr_inits2(long_prec, mp_myexp_long, mp_x, mp_exp_long, mp_myexp_relerr, mp_myexp_long_relerr, mp_myexp_rf_relerr, (mpfr_ptr)NULL);

	printf("x = "); while(scanf("%s", str_x) < 1);
	mpfr_set_str(mp_x, str_x, 10, MPFR_RNDN);

	mpfr_exp(mp_exp_long, mp_x, MPFR_RNDN);
/*	mpfr_myexp(mp_myexp, mp_x, n, MPFR_RNDN);

	mpfr_printf("myexp(%ld bits) = %RNe\n", prec, mp_exp);
	mpfr_printf("  exp(%ld bits) = %RNe\n", prec * 2, mp_exp_long);
	relerr(mp_exp_relerr, mp_exp, mp_exp_long);
	mpfr_printf("Relerr(pi) = %10.3RNe\n", mp_exp_relerr);
*/

	printf(" #term   Total     TErr      RErr\n");
	for(i = 1; i <= n; i++)
	{
		mpfr_myexp(mp_myexp, mp_x, i, MPFR_RNDN);
		mpfr_myexp(mp_myexp_long, mp_x, i, MPFR_RNDN);
		relerr(mp_myexp_relerr, mp_myexp, mp_exp_long);
		relerr(mp_myexp_long_relerr, mp_myexp_long, mp_exp_long);
		relerr(mp_myexp_rf_relerr, mp_myexp, mp_myexp_long);
		mpfr_printf("%5d %10.3RNe %10.3RNe %10.3RNe\n", i, mp_myexp_relerr, mp_myexp_long_relerr, mp_myexp_rf_relerr);
	}


	mpfr_clears(mp_myexp, mp_myexp_long, mp_exp_long, mp_myexp_relerr, mp_x, mp_myexp_long_relerr, mp_myexp_rf_relerr, (mpfr_ptr)NULL);

	return 0;
}
