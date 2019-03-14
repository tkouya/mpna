#include <stdio.h>

#include "mpfr.h"

// MPC library
#include "mpc.h"

int main()
{
	mpc_t a, b, c;
	mpfr_t real_num, imag_num;

	mpfr_set_default_prec(128);

	mpc_init2(a, mpfr_get_default_prec());
	mpc_init2(b, mpfr_get_default_prec());
	mpc_init2(c, mpfr_get_default_prec());

	mpfr_init(real_num); mpfr_init(imag_num);

	// a = sqrt(2) + sqrt(3) * i
	mpfr_sqrt_ui(real_num, 2UL, MPFR_RNDN);
	mpfr_sqrt_ui(imag_num, 3UL, MPFR_RNDN);
	mpc_set_fr_fr(a, real_num, imag_num, MPC_RNDNN);

	// b = -sqrt(5) + pi * i
	mpfr_sqrt_ui(real_num, 5UL, MPFR_RNDN);
	mpfr_neg(real_num, real_num, MPFR_RNDN);
	mpfr_const_pi(imag_num, MPFR_RNDN);
	mpc_set_fr_fr(b, real_num, imag_num, MPC_RNDNN);

	// print
	mpfr_printf("a = %40.32RNe + %40.32RNe * i\n", mpc_realref(a), mpc_imagref(a));
	mpfr_printf("b = %40.32RNe + %40.32RNe * i\n", mpc_realref(b), mpc_imagref(b));

	// Addition
	mpc_add(c, a, b, MPC_RNDNN);
	mpfr_printf("a + b = %40.32RNe + %40.32RNe * i\n", mpc_realref(c), mpc_imagref(c));

	// Subtraction
	mpc_sub(c, a, b, MPC_RNDNN);
	mpfr_printf("a - b = %40.32RNe + %40.32RNe * i\n", mpc_realref(c), mpc_imagref(c));

	// Multiplication
	mpc_mul(c, a, b, MPC_RNDNN);
	mpfr_printf("a * b = %40.32RNe + %40.32RNe * i\n", mpc_realref(c), mpc_imagref(c));

	// Division
	mpc_div(c, a, b, MPC_RNDNN);
	mpfr_printf("a / b = %40.32RNe + %40.32RNe * i\n", mpc_realref(c), mpc_imagref(c));

	// Sqrt
	mpc_sqrt(c, a, MPC_RNDNN);
	mpfr_printf("sqrt(a) = %40.32RNe + %40.32RNe * i\n", mpc_realref(c), mpc_imagref(c));

	// clear
	mpc_clear(a); mpc_clear(b); mpc_clear(c);
	mpfr_clear(real_num); mpfr_clear(imag_num);

	return 0;
}
