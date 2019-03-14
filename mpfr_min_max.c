#include <stdio.h>
#include <math.h>
#include <float.h>

#include "mpfr.h"

int main(void)
{
	long int prec;
	mpfr_t a, b, c, eps;

	mpfr_set_default_prec(128);

	mpfr_inits(a, b, c, eps, NULL);

	// PREC_MIN, MAX
	printf("MPFR_PREC_MIN, MAX = %ld, %ld\n", (long int)MPFR_PREC_MIN, (long int)MPFR_PREC_MAX);
	printf("exp_min, max = %ld, %ld\n", mpfr_get_emin(), mpfr_get_emax());

	// nan, +inf, -inf
	mpfr_set_nan(a);
	mpfr_set_inf(b, +1);
	mpfr_set_inf(c, -1);

	mpfr_printf("NAN, +INF, -INF = %RNe, %RNe, %RNe\n", a, b, c);

	// machine epsilon
	mpfr_set_ui(eps, 1UL, MPFR_RNDN);
	mpfr_nextabove(eps);
	mpfr_sub_ui(eps, eps, 1UL, MPFR_RNDN);

	mpfr_printf("prec = %ld, Machine epsilon = %RNe\n", mpfr_get_prec(eps), eps);

	// all precision
	printf("prec(in bits), machine epsilon\n");
	for(prec = MPFR_PREC_MIN; prec <= MPFR_PREC_MAX; prec *= 2)
	{
		// set precision
		mpfr_set_prec(eps, prec);

		// machine epsilon
		mpfr_set_ui(eps, 1UL, MPFR_RNDN);
		mpfr_nextabove(eps);
		mpfr_sub_ui(eps, eps, 1UL, MPFR_RNDN);

		if(mpfr_get_exp(eps) <= mpfr_get_emin())
			break;

		mpfr_printf("%15ld, %15.7RNe\n", mpfr_get_prec(eps), eps);
	}

	mpfr_clears(a, b, c, eps, NULL);

	return 0;
}
