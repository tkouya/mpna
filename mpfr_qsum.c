#include <stdio.h>
#include <stdlib.h>
#include "mpfr.h"


// quick-two-sum
// ret := a + b
// err := error of ret
void quick_two_sum(mpfr_t ret, mpfr_t err, mpfr_t a, mpfr_t b, mpfr_rnd_t rounding_mode)
{
	mpfr_add(ret, a, b, rounding_mode);
	mpfr_sub(err, ret, a, rounding_mode);
	mpfr_sub(err, b, err, rounding_mode);
}


int main()
{
	mpfr_rnd_t default_rounding_mode = MPFR_RNDN;
	mpfr_prec_t default_prec = 8;
	mpfr_t a, b, c, d;
	int base;

	mpfr_init2(a, default_prec); mpfr_set_ui(a, 2UL, default_rounding_mode); mpfr_sqrt(a, a, default_rounding_mode);
	mpfr_init2(b, default_prec); mpfr_set_ui(b, 3UL, default_rounding_mode); mpfr_sqrt(b, b, default_rounding_mode);
	mpfr_init2(c, default_prec);
	mpfr_init2(d, default_prec);

	// set a, b
	base = 2;
	printf("Input binary fp number: \n");
	printf("a = "); mpfr_inp_str(a, stdin, base, default_rounding_mode);
	mpfr_printf("a = %.7RNb\n", a); //mpfr_dump(a);

	printf("b = "); mpfr_inp_str(b, stdin, base, default_rounding_mode);
	mpfr_printf("b = %.7RNb\n", b);

	quick_two_sum(c, d, a, b, default_rounding_mode);

	mpfr_printf("a + b = %.7RNb\n", c);
	mpfr_printf("error = %.7RNb\n", d);

	mpfr_add(a, c, d, default_rounding_mode);

	mpfr_printf("ret + error = %.7RNb\n", a);

	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);
	mpfr_clear(d);

	return EXIT_SUCCESS;
}
