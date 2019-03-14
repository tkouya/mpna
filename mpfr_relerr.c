#include <stdio.h>

#include "mpfr.h"

void relerr(mpfr_t c, mpfr_t a, mpfr_t ad)
{
	mpfr_sub(c, a, ad, MPFR_RNDN);
	mpfr_div(c, c, ad, MPFR_RNDN);
	mpfr_abs(c, c, MPFR_RNDN);
}

int main()
{
	unsigned long prec;
	mpfr_t a, b, c;
	mpfr_t ad, bd, cd;

	printf("Input prec in bits: "); scanf("%ld", &prec);

	mpfr_set_default_prec(prec); // in bits

	// a = sqrt(2.0);
	// b = sqrt(3.0);
	mpfr_init_set_ui(a, 2UL, MPFR_RNDN); mpfr_sqrt(a, a, MPFR_RNDN);
	mpfr_init_set_ui(b, 3UL, MPFR_RNDN); mpfr_sqrt(b, b, MPFR_RNDN);
	mpfr_init(c);

	mpfr_printf("a = %50.43RNe\n", a);
	mpfr_printf("b = %50.43RNe\n", b);

	mpfr_init2(ad, prec * 2); mpfr_sqrt_ui(ad, 2UL, MPFR_RNDN);
	mpfr_init2(bd, prec * 2); mpfr_sqrt_ui(bd, 3UL, MPFR_RNDN);
	mpfr_init2(cd, prec * 2);

	relerr(c, a, ad);
	mpfr_printf("relerr(a) = %10.3RNe\n", c);
	relerr(c, b, bd);
	mpfr_printf("relerr(b) = %10.3RNe\n", c);

	mpfr_add(c, a, b, MPFR_RNDN);
	mpfr_add(cd, ad, bd, MPFR_RNDN);
	relerr(c, c, cd);
	mpfr_printf("relerr(a + b) = %10.3RNe\n", c);

	mpfr_sub(c, a, b, MPFR_RNDN);
	mpfr_sub(cd, ad, bd, MPFR_RNDN);
	relerr(c, c, cd);
	mpfr_printf("relerr(a - b) = %10.3RNe\n", c);

	mpfr_mul(c, a, b, MPFR_RNDN);
	mpfr_mul(cd, ad, bd, MPFR_RNDN);
	relerr(c, c, cd);
	mpfr_printf("relerr(a * b) = %10.3RNe\n", c);

	mpfr_div(c, a, b, MPFR_RNDN);
	mpfr_div(cd, ad, bd, MPFR_RNDN);
	relerr(c, c, cd);
	mpfr_printf("relerr(a * b) = %10.3RNe\n", c);

	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);

	mpfr_clear(ad);
	mpfr_clear(bd);
	mpfr_clear(cd);

	return 0;
}
