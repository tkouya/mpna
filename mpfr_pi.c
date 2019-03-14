#include <stdio.h>

#include "mpfr.h"

// Relative error
void relerr(mpfr_t c, mpfr_t a, mpfr_t ad)
{
	mpfr_sub(c, a, ad, MPFR_RNDN);
	mpfr_div(c, c, ad, MPFR_RNDN);
	mpfr_abs(c, c, MPFR_RNDN);
}

int main()
{
	mpfr_prec_t prec;
	mpfr_t mp_pi, mp_pi_long, mp_pi_relerr;

	printf("Input prec in bits: "); scanf("%ld", &prec);

	printf("prec = %ld\n", prec);
	mpfr_init2(mp_pi, prec);
	mpfr_init2(mp_pi_long, prec * 2);
	mpfr_init2(mp_pi_relerr, prec * 2);

	mpfr_const_pi(mp_pi, MPFR_RNDN);
	mpfr_const_pi(mp_pi_long, MPFR_RNDN);
	mpfr_printf("pi(%ld bits) = %RNe\n", prec, mp_pi);
	relerr(mp_pi_relerr, mp_pi, mp_pi_long);
	mpfr_printf("Relerr(pi) = %10.3RNe\n", mp_pi_relerr);

	mpfr_clear(mp_pi);
	mpfr_clear(mp_pi_long);
	mpfr_clear(mp_pi_relerr);

	return 0;
}
